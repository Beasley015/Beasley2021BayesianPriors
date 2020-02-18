######################################################
# Simulating data for augmentation & informed priors #
# Fall 2019                                          #
######################################################

# Setup -----------------------------------------------
# Load packages
library(R2jags)
library(tidyverse)
library(MASS)
library(abind)
library(boot)

# Set seed
set.seed(15)

# Global variables
nspec <- 15
nmiss <- 2 # Species present but not detected during sampling
naug <- 3 # Hypothetical species not exposed to sampling
ems <- nmiss + naug # Number of all-zero histories to add to observed data
nsite <- 30
nsurvey <- 4

Ks <- rep(nsurvey, nsite)

# Matrix of covariate responses
resp2cov <- rnorm(n = nspec+naug, sd = 0.25)

resp2cov <- sample(resp2cov)

# Covariate values for sites
cov <- sort(rnorm(n = nsite))

# Simulate occupancy data -------------------------------------
# Get probs from a beta distribution

sim.occ <- rbeta(n = nspec+nmiss, shape1 = 2, shape2 = 1)

# Write function to simulate true occupancy state
tru.mats <- function(spec=nspec+nmiss, site=nsite, alpha1=resp2cov){
  #Get site-level psi to account for covariates
  alpha0 <- logit(sim.occ)
  
  logit.psi <- matrix(NA, nrow = spec, ncol = site)
  
  for(i in 1:spec){
    logit.psi[i,] <- alpha0[i] + alpha1[i]*cov
  }
  
  psi <- plogis(logit.psi)
  
  #create list of abundance vectors
  nlist<-list()
  for(a in 1:spec){
    nlist[[a]] <- rbinom(n = site, size = 1, prob = psi[a,])
  }
  
  #Turn abundance vectors into abundance matrix
  ns<-do.call(rbind, nlist)
  
  return(list(ns, psi))
}

tru <- tru.mats()[[1]]
mu.psi <- rowMeans(tru.mats()[[2]])

# Simulate detection process ----------------------------------
# Load detections from Master's work
masters.mod <- readRDS("modelsampledglades.rds")
masters.v <- masters.mod$sims.list$v

spec.p <- plogis(masters.v[,1:8])

spec.det <- colMeans(spec.p)

# Get maximum likelihood estimates for params of beta distribution
bet.det <- fitdistr(x = spec.det, start = list(shape1 = 1, shape2 = 1), "beta")

# Reorder true occurrence matrix, psi, and cov responses
reorder <- function(x){
  if(length(dim(x)) == 0){
    pos <- c(which(psi == min(psi)), which(psi == median(psi)))
    nondet <- x[pos]
    x <- x[-pos]
    x <- c(x, nondet)
  } else {
    pos <- c(which(psi == min(psi)), which(psi == median(psi)))
    nondet <- x[pos,]
    x <- x[-pos,]
    x <- rbind(x, nondet)
  }
  return(x)
}

tru <- reorder(x = tru)
resp2cov <- reorder(x = resp2cov)
occ <- reorder(x = sim.occ)

# Generate detection probabilities from beta dist with above params
sim.dets <- rbeta(n = nspec, shape1 = bet.det$estimate[1], shape2 = bet.det$estimate[2])

# Assign selected species a detection probability of 0
sim.dets <- c(sim.dets, rep(0,2))

# Function to create encounter histories
trap.hist <- function(mat, det, specs=nspec, sites=nsite, survs=nsurvey){
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

obs.data <- trap.hist(mat = tru, det = sim.dets)

# Augment the observed dataset ------------------------------------
ems.array <- array(0, dim = c(nsite, nsurvey, ems))
obs.aug <- abind(obs.data, ems.array, along = 3)

# Add prior information --------------------------------
uninf <- list(rep(0, nspec+ems), rep(0, nspec+ems))
weakinf <- list(c(rep(0, nspec), logit(sim.occ[16:17])*0.1, rep(0, naug)),
                  c(rep(0, nspec), round(resp2cov[16:17])*0.1, rep(0, naug)))
modinf <- list(c(rep(0, nspec), logit(sim.occ[16:17])*0.5, rep(0, naug)),
               c(rep(0, nspec), round(resp2cov[16:17])*0.5, rep(0, naug)))
weakmisinf <- list(c(rep(0, nspec), logit(sim.occ[16:17])*-0.1, rep(0, naug)),
                   c(rep(0, nspec), round(resp2cov[16:17])*-0.1, rep(0, naug)))
modmisinf <- list(c(rep(0, nspec), logit(sim.occ[16:17])*-0.5, rep(0, naug)),
                  c(rep(0, nspec), round(resp2cov[16:17])*-0.5, rep(0, naug)))

# Write Models ----------------------------
# Model without augmentation
cat("
    model{
      
    # Define hyperprior distributions: intercepts

    #Intercepts
    mean.a0 ~ dunif(0,1)
    a0.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mean.a1 ~ dunif(0,1)
    a1.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a1 ~ dgamma(0.1, 0.1)
    
    mean.b0 ~ dunif(0,1)
    b0.mean <- log(mean.b0)-log(1-mean.b0)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    for(i in 1:spec){
    #create priors from distributions above
    
    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i]
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
    
    N<-spec
    
    }
    ", file = "noaug.txt")

# Model with augmentation
cat("
    model{
      
    # Define hyperprior distributions: intercepts
    omega ~ dunif(0,1)
    
    #Intercepts
    mean.a0 ~ dunif(0,1)
    a0.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mean.a1 ~ dunif(0,1)
    a1.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a1 ~ dgamma(0.1, 0.1)
    
    mean.b0 ~ dunif(0,1)
    b0.mean <- log(mean.b0)-log(1-mean.b0)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    for(i in 1:(spec+aug)){
    # Create priors from hyperpriors
    w[i] ~ dbern(omega)
    #indicates whether or not species is exposed to sampling
    
    a0[i] ~ dnorm(a0.mean+info1[i], tau.a0)
    a1[i] ~ dnorm(a1.mean+info2[i], tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i] * w[i]
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
    n0<-sum(w[(spec+1):(spec+aug)])
    N<-spec+n0
    
    }
    ", file = "aug_model.txt")

# Write function for sending model to gibbs sampler --------------------------------
VivaLaMSOM <- function(J, K, obs, spec, aug = 0, cov, textdoc, info1 = NULL, 
                       info2 = NULL, burn = 2500, iter = 8000, thin = 10){
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, cov = cov)
  if(textdoc == 'aug_model.txt'){
    datalist$aug <- aug
  }
  
  if(is.null(info1) == F){
    datalist$info1 <- info1
    datalist$info2 <- info2
  }

  # Specify parameters
  parms <- c('N', 'a0.mean', 'b0.mean', 'a0', 'b0', 'a1', 'Z')

  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
    omega.guess <- runif(1,0,1)
    mu.psi.guess <- runif(1, 0.25, 1)
    inits <- list(
         a0 = rnorm(n = (spec+aug)), a1 = rnorm(n = (spec+aug)),
         b0 = rnorm(n = (spec+aug)),
         Z = maxobs
    )
    if(textdoc == 'aug_model.txt'){
      inits$omega <- omega.guess
      inits$w <- c(rep(1,spec), rbinom(n = aug, size=1, prob=omega.guess))
    }
    
    return(inits)
  }

  #JAGS command
  model <- jags(model.file = textdoc, data = datalist, n.chains = 3,
                parameters.to.save = parms, inits = init.values, n.burnin = burn,
                n.iter = iter, n.thin = thin)
    
    return(model)
}

# Run sims ------------------------------------
# mod.noaug <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.data, cov = cov, spec = nspec,
#            textdoc = 'noaug.txt')
# saveRDS(mod.noaug, file = "mod_noaug.rds")

# mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec, 
#            textdoc = 'aug_model.txt', aug = nmiss+naug, info1 = uninf[[1]],
#            info2 = uninf[[2]], burn = 2500, iter = 10000, thin = 10)
# saveRDS(mod.uninf, file = "mod_uninf.rds")

# mod.inf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec, 
#                       textdoc = 'aug_model.txt', aug = nmiss+naug, 
#                       info1 = weakinf[[1]], info2 = weakinf[[2]], burn = 5000, 
#                       iter = 12000, thin = 5)
# saveRDS(mod.inf.weak, file = "mod_inf_weak.rds")

# mod.inf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec, 
#                       textdoc = 'aug_model.txt', aug = nmiss+naug, 
#                       info1 = modinf[[1]], info2 = modinf[[2]], burn = 7000, 
#                       iter = 12000, thin = 3)
# saveRDS(mod.inf, file = "mod_inf.rds")

mod.misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, 
                              spec = nspec, textdoc = 'aug_model.txt', 
                              aug = nmiss+naug, info1 = weakmisinf[[1]], 
                              info2 = weakmisinf[[2]], burn = 5000, iter = 10000, 
                              thin = 5)

mod.misinf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
                         textdoc = 'aug_model.txt', aug = nmiss+naug, 
                         info1 = modmisinf[[1]], info2 = modmisinf[[2]], burn = 2500,
                         iter = 10000, thin = 10)

# Function to compare estimates of N -----------------------------------
getmode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

Ns <- as.vector(mod.uninf$BUGSoutput$sims.list$N)
Ns %>%
  table() %>%
  data.frame() %>%
  {. ->> ns.frame}
colnames(ns.frame) <- c("N_Species", "Freq")

Ns.mode <-getmode(Ns)
Ns.mean <- mean(Ns)

ggplot(data = ns.frame, aes(x = as.integer(as.character(N_Species)), y = Freq))+
  geom_col()+
  geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"), size = 1.5)+
  geom_vline(aes(xintercept = Ns.mean, linetype = "Estimated"), size = 1.5)+
  scale_linetype_discrete(name = "")+
  labs(x = "Number of Species")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic(base_size = 18)+
  theme(axis.text.y = element_blank(), legend.key.size = unit(2.5, "lines"))
