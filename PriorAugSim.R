######################################################
# Simulating data for augmentation & informed priors #
# Fall 2019                                          #
######################################################

# Setup -----------------------------------------------
# Load packages
library(R2OpenBUGS)
library(vcdExtra)
library(tidyverse)
library(MASS)

# Set seed
set.seed(15)

# Global variables
nspec <- 15
naug <- 1
nsite <- 30
nsurvey <- 4

Ks <- rep(nsurvey, nsite)

# Matrix of covariate responses
resp2cov <- c(rnorm(n = 9, mean = 1, sd = 0.25), 
                 rnorm(n = 4, mean = 0, sd = 0.25),
                 rnorm(n = 3, mean = -1, sd = 0.25))

# Covariate values for sites
cov.bin <- as.numeric(rbernoulli(n = nsite, 0.5))
cov.cont <- rnorm(n = nsite, mean = 0, sd = 2)

cov.vals <- cbind(cov.bin, cov.cont)
cov.names <- c("bin", "cont")

# Function for abundance data -------------------------------------
tru.mats <- function(spec = nspec + naug, site = nsite, cov, resp){
  #Draw lambdas from a logseries distribution
  mean.lambdas <- rgamma(spec, shape=1.5, rate=2)
  
  #Get site-level lambdas to account for covariates
  alpha0 <- log(mean.lambdas) #log-scale intercept
  alpha1 <- exp(resp) #log-scale covariate responses
  
  log.lambdas <- matrix(NA, nrow = spec, ncol = site)
  
  for(i in 1:spec){
    log.lambdas[i,] <- alpha0[i] + alpha1[i]*cov
  }

  lambdas <- exp(log.lambdas)  #inverse link transformation
  
  #create list of abundance vectors
  nlist<-list()
  for(a in 1:spec){
    nlist[[a]] <- rpois(n = site, lambda = lambdas[a])
  }
  
  #Turn abundance vectors into abundance matrix
  ns<-do.call(rbind, nlist)
  
  ns[ns > 1] <- 1
  
  output <- list(ns, mean.lambdas)
  return(output)
}

#Get true abundances
for(i in 1:ncol(cov.vals)){
    assign(paste(cov.names[i], sep = ""),
          tru.mats(cov = cov.vals[,i], resp = resp2cov))
}

trus <- list(bin[[1]], cont[[1]])
lambdas <- list(bin[[2]], cont[[2]])

#name objects in objects to keep them straight
sim.names <- c('bin', 'cont')

names(trus) <- sim.names
names(lambdas) <- sim.names

# Simulate detection process ----------------------------------
# Load model results from Master's work
masters.mod <- readRDS("modelsampledglades.rds")

masters.v <- masters.mod$sims.list$v

spec.p <- plogis(masters.v[,1:8])

spec.det <- colMeans(spec.p)

# Get maximum likelihood estimates for params of beta distribution
bet <- fitdistr(x = spec.det, start = list(shape1 = 1, shape2 = 1), "beta")

# Rare species is undetected
for(i in 1:length(trus)){
  row <- trus[[i]][which(lambdas[[i]] == min(lambdas[[i]])),]
  trus[[i]] <- trus[[i]][-which(lambdas[[i]] == min(lambdas[[i]])),]
  trus[[i]] <- rbind(trus[[i]], row)
}

# DO THE SAME FOR COVARIATE VALUES
reorder.covs <- function(lbd, FUN){
    resp <- resp2cov[which(lbd == FUN(lbd))]
    resp2cov <- resp2cov[-which(lbd == FUN(lbd))]
    resp2cov <- c(resp2cov, resp)
}

resp.bin <- reorder.covs(lambdas$bin, FUN = min)
resp.cont <- reorder.covs(lambdas$cont, FUN = min)

# Common species is undetected
# for(i in 1:length(trus)){
#   row <- trus[[i]][which(lambdas[[i]] == median(lambdas[[i]])),]
#   trus[[i]] <- trus[[i]][-which(lambdas[[i]] == median(lambdas[[i]])),]
#   trus[[i]] <- rbind(trus[[i]], row)
# }

# resp.bin <- reorder.covs(lambdas$bin, FUN = median)
# resp.cont <- reorder.covs(lambdas$cont, FUN = median)

# Generate detection probabilities from beta dist with above params
sim.dets <- rbeta(n = nspec, shape1 = bet$estimate[1], shape2 = bet$estimate[2])

# Assign last species a detection probability of 0
sim.dets <- c(sim.dets, 0)

# Function to create encounter histories
trap.hist <- function(mat, det, specs=nspec+naug, sites=nsite, survs=nsurvey){
  #Detection intercept and cov responses
  beta0<-qlogis(det) #put it on logit scale
  #Responses to a cov would go here
  
  #Logit link function
  logit.p <- array(NA, dim = c(sites, survs, specs))
  for(i in 1:specs){
    logit.p[,,i] <- beta0[i]
  }
  
  p <- plogis(logit.p)
  
  #Simulate observation data
  L<-list()
  
  for(b in 1:specs){
    y<-matrix(NA, ncol = survs, nrow = sites)
    for(a in 1:survs){
      y[,a]<-rbinom(n = sites, size = mat[b,], prob = p[,,b])
    }
    L[[b]]<-y
  }
  
  #Smash it into array
  obsdata<-array(as.numeric(unlist(L)), dim=c(sites, survs, specs))
  
  return(obsdata)
}

obs.data <- lapply(trus, trap.hist, det = sim.dets)

names(obs.data) <- sim.names

# Write Models ----------------------------
# Model without augmentation
cat("
    model{
      
    # Define hyperprior distributions: intercepts
    
    omega ~ dunif(0,1)
    
    #Intercepts
    a0.mean ~ dnorm(0,0.001)
    sigma.a0 ~ dunif(0,10)
    tau.a0 <- 1/(sigma.a0*sigma.a0)
    
    a1.mean ~ dnorm(0,0.001)
    sigma.a1 ~ dunif(0,10)
    tau.a1 <- 1/(sigma.a1*sigma.a1)
    
    b0.mean ~ dnorm(0,0.001)
    sigma.b0 ~ dunif(0,10)
    tau.b0 <- 1/(sigma.b0*sigma.b0)
    
    for(i in 1:nspec){
    #create priors from distributions above
    w[i] ~ dbern(omega)
    #indicates whether or not species is exposed to sampling
    
    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i]*w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    #Estimate detection of i at point j during sampling period k
    for(k in 1:K[j]){
    logit(p[j,k,i]) <-  b0[i]
    mu.p[j,k,i] <- p[j,k,i]*Z[j,i] 
    #The addition of Z means that detecting a species depends on its occupancy
    obs[j,k,i] ~ dbern(mu.p[j,k,i])
    }
    }
    }
    
    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
    n0<-sum(w[(nspec+1):(nspec+naug)])
    N<-nspec+n0
    
    }
    ", file = "noaug.txt")

# Model with augmentation, uninformed priors
cat("
    model{
      # Define hyperprior distributions: intercepts
    
    omega ~ dunif(0,1)
    
    #Intercepts
    a0.mean ~ dnorm(0,0.001)
    sigma.a0 ~ dunif(0,10)
    tau.a0 <- 1/(sigma.a0*sigma.a0)

    a1.mean ~ dnorm(0,0.001)
    sigma.a1 ~ dunif(0,10)
    tau.a1 <- 1/(sigma.a1*sigma.a1)
    
    b0.mean ~ dnorm(0,0.001)
    sigma.b0 ~ dunif(0,10)
    tau.b0 <- 1/(sigma.b0*sigma.b0)
    
    for(i in 1:nspec+naug){
    #create priors from distributions above
    w[i] ~ dbern(omega)
    #indicates whether or not species is exposed to sampling
    
    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i]*w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    #Estimate detection of i at point j during sampling period k
    for(k in 1:K[j]){
    logit(p[j,k,i]) <-  b0[i]
    mu.p[j,k,i] <- p[j,k,i]*Z[j,i] 
    #The addition of Z means that detecting a species depends on its occupancy
    obs[j,k,i] ~ dbern(mu.p[j,k,i])
    }
    }
    }
    
    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
    n0<-sum(w[(nspec+1):(nspec+naug)])
    N<-nspec+n0
    
    }
    ", file = "aug_uninf.txt")

# Model with informed priors
cat("
    ")

# Model with misinformed priors
cat("
    ")

# Write function for sending model to gibbs sampler --------------------------------
VivaLaMSOM <- function(J, K, obs, nspec, naug=NULL, textdoc, covar){
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, nspec = nspec, naug = naug)
  if(textdoc == "cov.txt"){
    list.append(datalist, cov=covar)
  }

  # Specify parameters
  parms <- list('Z', 'N', 'a0', 'b0', 'mu.psi', 'mu.p')
  if(textdoc == "cov.txt"){
   list.append(parms, 'a1') 
  }

  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
    omega.guess <- runif(1,0,1)
    list(omega = omega.guess,
         w=c(rep(1,nspec), rbinom(n = naug, size=1, prob=omega.guess)),
         a0 = rnorm(n = (nspec+naug)), a1 = rnorm(n = (nspec+naug)),
         b0 = rnorm(n = (nspec+naug)),
         Z = maxobs)
  }

  # JAGS command: this isn't working; hates my initial values for some reason
  # model <- jags(model.file = "nocov.txt", data = datalist, n.chains = 3,
  #               parameters.to.save = parms, inits = init.values, n.burnin = 100,
  #               n.iter = 1000)

  # BUGS works though
    model <- bugs(model.file = textdoc, data = datalist, n.chains = 3,
                  parameters.to.save = parms, inits = init.values, n.burnin = 10,
                  n.iter = 100, debug = F)
    
    return(model)
}

# Run sims -------------------------------------
outputs.nocov <- list()
for(i in 1:length(obs.data)){
  outputs.nocov[[i]] <- VivaLaMSOM(J=nsite, K=Ks, nspec=nspec, naug=naug, 
                             obs = obs.data[[i]], textdoc = "nocov.txt")  
  print(i)
}
names(outputs.nocov) <- sim.names

outputs.bin.uninf <- list()

outputs.cont.uninf <- list()

outputs.bin.inf <- list()

outputs.cont.inf <- list()
