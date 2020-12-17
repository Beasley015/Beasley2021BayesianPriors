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
library(viridis)
library(patchwork)
library(fitdistrplus)
library(grid)

# Set seed
set.seed(16)

# Global variables
nspec <- 20
nmiss <- 2 # Species present but not detected during sampling
nsite <- 30
nsurvey <- 4

Ks <- rep(nsurvey, nsite)

# Vector of covariate responses
resp2cov <- c(rnorm(n = 7, sd = 0.25),
              rnorm(n = 8, mean = 3, sd = 0.25),
              rnorm(n = 7, mean = -3, sd = 0.25))

resp2cov <- sample(resp2cov)

# Covariate values for sites
cov <- sort(rnorm(n = nsite))

# Simulate occupancy data -------------------------------------
# Get probs from a beta distribution
sim.occ <- rbeta(n = nspec+nmiss, shape1 = 2, shape2 = 3)

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
psi <- rowMeans(tru.mats()[[2]])

# Simulate detection process ---------------------------------
# Generate mean detection probabilities from beta dist
mean.p <- rbeta(n = nspec+nmiss, shape1 = 2, shape2 = 8)
mean.p <- sort(mean.p, decreasing = T)

# Generate detection histories
get.obs <- function(mat, specs){
  #Detection intercept and cov responses
  beta0<-logit(mean.p) #put it on logit scale

  #Logit link function
  logit.p <- array(NA, dim = c(nsite, nsurvey, specs))
  for(i in 1:specs){
    for(j in 1:nsite){
      for(k in 1:nsurvey){
        logit.p[j,,i] <- beta0[i]
      }
    }
  }

  p <- plogis(logit.p)

  #Simulate observation data
  L<-list()

  for(b in 1:specs){
    y<-matrix(NA, ncol = nsite, nrow = nsurvey)
    for(a in 1:nsurvey){
      y[a,]<-rbinom(n = nsite, size = 1, prob = p[,,b]*mat[b,])
    }
    L[[b]]<-t(y)
  }

  #Smash it into array
  obsdata<-array(as.numeric(unlist(L)), dim=c(nsite, nsurvey, specs))

  #Get site-level p
  ps <- apply(p, c(1,3), mean)
  
  #List of outputs
  out.list <- list(obs.data = obsdata, ps = ps)

  return(out.list)
}

obs.data <- get.obs(mat = tru, specs = nspec+nmiss)[[1]]
p <- get.obs(mat = tru, specs = nspec+nmiss)[[2]]

#Checkpoint: get observed data matrix
maxobs <- apply(obs.data, c(1,3), max)
colSums(maxobs)

# Remove species with fewest detections: these will be "undetected" species
obs.data <- obs.data[,,-c(21:22)]

# Function to reorder true values, if needed
# reorder <- function(x){
#   if (length(dim(x)) == 0){
#     nondets <- which(colSums(maxobs) == 0)
#     copy <- x[nondets]
#     x <- x[-nondets]
#     new <- c(x, copy)
#     return(new)
#     }
#   else {
#     nondets <- which(colSums(maxobs) == 0)
#     copy <- x[nondets,]
#     x <- x[-nondets,]
#     new <- rbind(x, copy)
#     return(new)
#     }
# }
# 
# sim.occ <- reorder(sim.occ)
# mean.p <- reorder(mean.p)
# 
# resp2cov <- reorder(resp2cov)
# 
# tru <- reorder(tru)

# Augment the observed dataset ------------------------------------
ems.array <- array(0, dim = c(nsite, nsurvey, nmiss))
obs.aug <- abind(obs.data, ems.array, along = 3)

# Add prior information --------------------------------
uninf <- "for(i in 1:(spec+aug)){
          #Create priors from hyperpriors
            w[i] ~ dbern(omega)
            
            a0[i] ~ dnorm(a0.mean, tau.a0)
            a1[i] ~ dnorm(a1.mean, tau.a1)

            b0[i] ~ dnorm(b0.mean, tau.b0)"


weakinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -4,0, 0)
            inf.mean1 <- c(0, 0,3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.85, 0.15)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,2] <- weights[2]/inf.var0[g[i]+1]
              lb0[i,1] <- weights[1]/(1/tau.a0)
              lb1[i,1] <- weights[1]/(1/tau.a1)
              lb1[i,2] <- weights[2]/inf.var1[g[i]+1]
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, (1/pooled.var0[i]), 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, (1/pooled.var1[i]), 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

modinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -4,0, 0)
            inf.mean1 <- c(0, 0,3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.5, 0.5)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,2] <- weights[2]/inf.var0[g[i]+1]
              lb0[i,1] <- weights[1]/(1/tau.a0)
              lb1[i,1] <- weights[1]/(1/tau.a1)
              lb1[i,2] <- weights[2]/inf.var1[g[1]+1]
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, 1/pooled.var0[i], 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, 1/pooled.var1[i], 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

stronginf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -4,0, 0)
            inf.mean1 <- c(0, 0,3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.15, 0.85)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,2] <- weights[2]/inf.var0[g[i]+1]
              lb0[i,1] <- weights[1]/(1/tau.a0)
              lb1[i,1] <- weights[1]/(1/tau.a1)
              lb1[i,2] <- weights[2]/inf.var1[g[i]+1]
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, 1/pooled.var0[i], 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, 1/pooled.var1[i], 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"
  
weakmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, 3,3, 0)
            inf.mean1 <- c(0, -3,-3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.85, 0.15)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,2] <- weights[2]/inf.var0[g[i]+1]
              lb0[i,1] <- weights[1]/(1/tau.a0)
              lb1[i,1] <- weights[2]/inf.var1[g[i]+1]
              lb1[i,2] <- weights[1]/(1/tau.a1)
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, 1/pooled.var0[i], 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, 1/pooled.var1[i], 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"
 
modmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, 3,3, 0)
            inf.mean1 <- c(0, -3,-3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.5, 0.5)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,2] <- weights[2]/(1/tau.a0)
              lb0[i,1] <- weights[1]/inf.var0[g[i]+1]
              lb1[i,1] <- weights[1]/inf.var1[g[i]+1]
              lb1[i,2] <- weights[2]/(1/tau.a1)
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, 1/pooled.var0[i], 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, 1/pooled.var1[i], 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

strongmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, 3,3, 0)
            inf.mean1 <- c(0, -3,-3, 0)
            
            inf.var0 <- c(1, 0.5,0.5, 1)
            inf.var1 <- c(1, 0.5,0.5, 1)
            
            weights <- c(0.15, 0.85)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              w[i] ~ dbern(omega)
              
              lb0[i,1] <- weights[1]/(1/tau.a0)
              lb0[i,2] <- weights[2]/inf.var0[g[i]+1]
              lb1[i,2] <- weights[2]/inf.var1[g[i]+1]
              lb1[i,1] <- weights[1]/(1/tau.a1)
              
              pooled.var0[i] <- 1/sum(lb0[i,])
              pooled.var1[i] <- 1/sum(lb1[i,])
              
              pooled.mean0[i] <- sum(lb0[i,]*c(a0.mean,inf.mean0[g[i]+1]))
                                 *pooled.var0[i]
              pooled.mean1[i] <- sum(lb1[i,]*c(a1.mean,inf.mean1[g[i]+1]))
                                 *pooled.var1[i]
              
              a0[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean0[i], 
                            a0.mean), 
                            ifelse(i==21 || i==22, 1/pooled.var0[i], 
                            tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==21 || i==22, pooled.mean1[i], 
                            a1.mean), 
                            ifelse(i==21 || i ==22, 1/pooled.var1[i], 
                            tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

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
write.model <- function(priors){
  mod <- paste("
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

    ",priors,"

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
    ")
 writeLines(mod, "aug_model.txt") 
}

# Write function for sending model to gibbs sampler --------------------------------
VivaLaMSOM <- function(J, K, obs, spec, aug = 0, cov, textdoc, 
                       priors = uninf, burn = 2500, iter = 8000, 
                       thin = 10){
  # Write model for augmented datasets
  if(textdoc == 'aug_model.txt')
    write.model(priors = priors)
  
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, cov = cov)
  if(textdoc == 'aug_model.txt'){
    datalist$aug <- aug
  }

  # Specify parameters
  parms <- c('a0.mean', 'tau.a0', 'a1.mean', 'tau.a1', 'a0', 'psi',
             'b0', 'a1','Z', 'pooled.mean0', 'pooled.var0', 
             'pooled.mean1', 'pooled.var1', 'N')

  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
    inits <- list(
         a0 = rnorm(n = (spec+aug)),
         a1 = rnorm(n = (spec+aug)),
         b0 = rnorm(n = (spec+aug)),
         Z = maxobs
    )
    if(textdoc == 'aug_model.txt'){
      inits$w <- c(rep(1,spec), rbinom(n = aug, size=1, 
                                       prob = runif(1,0,1)))
    }
    
    return(inits)
  }

  #JAGS command
  model <- jags(model.file = textdoc, data = datalist, n.chains = 3,
                parameters.to.save = parms, inits = init.values, 
                n.burnin = burn, n.iter = iter, n.thin = thin)
    
  return(model)
}

# Run sims ------------------------------------
# mod.noaug <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.data, cov = cov,
#                         spec = nspec, textdoc = 'noaug.txt')
# saveRDS(mod.noaug, file = "mod_noaug.rds")
# 
# mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                         spec = nspec,textdoc = 'aug_model.txt',
#                         aug = nmiss, burn = 2500, iter = 10000,
#                         thin = 10)
# saveRDS(mod.uninf, file = "mod_uninf.rds")
# 
# inf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                            spec = nspec,textdoc = 'aug_model.txt',
#                            aug = nmiss, priors = weakinf, burn = 2500,
#                            iter = 10000, thin = 5)
# saveRDS(inf.weak, file = "inf_weak.rds")
# 
# inf.mod <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                       spec = nspec, textdoc = 'aug_model.txt',
#                       aug = nmiss, priors = modinf, burn = 7000,
#                       iter = 12000, thin = 3)
# saveRDS(inf.mod, file = "inf_mod.rds")
# 
# inf.strong <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                          spec = nspec, textdoc = 'aug_model.txt',
#                          aug = nmiss, priors = stronginf, burn = 7000,
#                          iter = 12000, thin = 3)
# saveRDS(inf.strong, file = "inf_strong.rds")
# 
# misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                           spec = nspec, textdoc = 'aug_model.txt',
#                           aug = nmiss, priors = weakmisinf, burn = 5000,
#                           iter = 10000, thin = 5)
# saveRDS(misinf.weak, file = "misinf_weak.rds")
# 
# misinf.mod <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                          spec = nspec, textdoc = 'aug_model.txt',
#                          aug = nmiss, priors = modmisinf, burn = 2500,
#                          iter = 10000, thin = 10)
# saveRDS(misinf.mod, file = "misinf_mod.rds")
# 
# misinf.strong <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                             spec = nspec, textdoc = 'aug_model.txt',
#                             aug = nmiss, priors = strongmisinf,
#                             burn = 2000, iter = 10000, thin = 5)
# saveRDS(misinf.strong, file = "misinf_strong.rds")

# Load models -------------------------------
noaug <- readRDS(file = "mod_noaug.rds")

# Augmented models
uninf <- readRDS("mod_uninf.rds")
inf.weak <- readRDS("inf_weak.rds")
inf.mod <- readRDS("inf_mod.rds")
inf.strong <- readRDS("inf_strong.rds")
misinf.weak <- readRDS("misinf_weak.rds")
misinf.mod <- readRDS("misinf_mod.rds")
misinf.strong <- readRDS("misinf_strong.rds")

# Put models in a list
mod.outputs <- list(inf.weak, inf.mod, inf.strong, misinf.weak, 
                    misinf.mod, misinf.strong)
names(mod.outputs) <- c("inf.weak", "inf.mod", "inf.strong", 
                        "misinf.weak", "misinf.mod", "misinf.strong")

# Compare priors and posteriors ------------
# Write function to pull priors/posteriors and plot
prior.agg0 <- function(x, spec.ID, inf.means){
  pooled.mean <- mean(x$BUGSoutput$sims.list$pooled.mean0[,spec.ID])
  pooled.sd <- mean(sqrt(x$BUGSoutput$sims.list$pooled.var0[,spec.ID]))

  inf.mean <- inf.means[spec.ID-20]
  inf.var <- 0.5
  
  comm.mean <- mean(x$BUGSoutput$sims.list$a0.mean)
  comm.sd <- mean(sqrt(1/x$BUGSoutput$sims.list$tau.a0))
  
  post.mean <- mean(x$BUGSoutput$sims.list$a0[,spec.ID])
  post.sd <- sd(x$BUGSoutput$sims.list$a0[,spec.ID])

  plot <- ggplot()+
    stat_function(fun = dnorm, n = 1000, 
                  args = list(mean = pooled.mean, sd = pooled.sd),
                  size = 1)+
    stat_function(fun = dnorm, n = 1000, 
                  args = list(mean = inf.mean, sd = sqrt(inf.var)),
                  size = 1, linetype = "dashed")+
    stat_function(fun = dnorm, n = 1000, 
                  args = list(mean = comm.mean, sd = comm.sd),
                  size = 1, linetype = "dotted")+
    stat_function(fun = dnorm, n = 1000, 
                  args = list(mean = post.mean, sd = post.sd),
                  size = 1, color = "red")+
    xlim(c(-6, 5))+
    labs(y = "Density")+
    theme_bw(base_size = 16)+
    theme(panel.grid = element_blank(), 
          axis.title.x = element_blank())
  
  return(plot)
}

# Create subplots for patchwork
a0.21.inf <- lapply(mod.outputs[1:3], prior.agg0, 
                    spec.ID = 21, inf.means = c(-4, 0))
a0.21.misinf <- lapply(mod.outputs[4:6], prior.agg0, spec.ID = 21, 
                       inf.means = c(3,3))

# Create patchwork objects
a021inf <- (a0.21.inf[[1]]+a0.21.inf[[2]]+a0.21.inf[[3]])
a021misinf <- (a0.21.misinf[[1]]+a0.21.misinf[[2]]+
                 a0.21.misinf[[3]])

# Create layout and add tags
a021patch <- a021inf/a021misinf

a021patch[[1]] <- a021patch[[1]]+plot_layout(tag_level = 'new') 
a021patch[[2]] <- a021patch[[2]]+plot_layout(tag_level = 'new')

# Final plot
a021plot <- a021patch+
  plot_annotation(tag_levels = c('A', '1'))

# Save that sonofabitch
# ggsave(filename = "a021plot.jpeg", a021plot, width = 8, 
#        height = 5)

# Do it again with the other species
# a0.22.inf <- lapply(mod.outputs[1:3], prior.agg0, 
#                     spec.ID = 22, inf.means = c(-4, 0))
# a0.22.misinf <- lapply(mod.outputs[4:6], prior.agg0, spec.ID = 22, 
#                        inf.means = c(3,3))
# 
# a022inf <- (a0.22.inf[[1]]+a0.22.inf[[2]]+a0.22.inf[[3]])
# a022misinf <- (a0.22.misinf[[1]]+a0.22.misinf[[2]]+
#                  a0.22.misinf[[3]])
# 
# a022patch <- a022inf/a022misinf
# 
# a022patch[[1]] <- a022patch[[1]]+plot_layout(tag_level = 'new') 
# a022patch[[2]] <- a022patch[[2]]+plot_layout(tag_level = 'new')
# 
# a022plot <- a022patch+
#   plot_annotation(tag_levels = c('A', '1'))

# Compare Ns -----------------
# Create list that includes uninformed priors
biglist <- list(uninf, inf.weak, inf.mod, 
                inf.strong, misinf.weak, misinf.mod, 
                misinf.strong)

# Write function to get Ns
get.ns <- function(jag){
  getmode <- function(x) {
    uniqx <- unique(x)
    uniqx[which.max(tabulate(match(x, uniqx)))]
  }
  
  Ns <- as.vector(jag$BUGSoutput$sims.list$N)
  Ns %>%
    table() %>%
    data.frame() %>%
    {. ->> ns.frame}
  colnames(ns.frame) <- c("N_Species", "Freq")
  
  Ns.mode <-getmode(Ns)
  Ns.mean <- mean(Ns)
  Ns.median <- median(Ns)
  
  Ns.plot <- ggplot(data = ns.frame, 
                    aes(x = as.integer(as.character(N_Species)),
                                       y = Freq))+
    geom_col(width = 0.95, color = 'lightgray')+
    geom_vline(aes(xintercept = Ns.median, linetype = "Estimated"),
               size = 1.5)+
    geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"),
               size = 1.5)+
    scale_linetype_manual(values = c("Estimated"="dotted", 
                                     "True"="solid"),
                          name = "", 
                          labels = c("Median Estimate", "True"))+
    labs(x = "Estimated Species", y = "Frequency")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic(base_size = 18)+
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(), 
          axis.title.x = element_blank(),
          legend.key.height = unit(40, units = 'pt'),
          aspect.ratio = 1/1)
  
  out.list <- list(plot = Ns.plot, mode = Ns.mode, mean = Ns.mean,
                   median = Ns.median)
  
  return(out.list)
}

N.outs <- lapply(biglist, get.ns)

# Put histograms in single figure
histos <- map(N.outs, 1) 

# Create patchwork objects
histos.uninf <- histos[[1]]+
  plot_annotation(tag_levels = "A")&
  theme(plot.margin = unit(c(5,5, 0, 5.5, 5.5), units = "point"))
histos.inf <- (histos[[2]]/histos[[3]]/histos[[4]])&
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), units = "point"))
histos.misinf <- (histos[[5]]/histos[[6]]/histos[[7]])&
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 0), units = "point"))

# Final plot
Ns.megaplot <- histos.uninf/(histos.inf|histos.misinf)+
  plot_annotation(tag_levels = "A")+
  plot_layout(guides = "collect", heights = c(1,3.5)) &
  xlim(19.5, 22.5)
  
# ggsave(Ns.megaplot, filename = "ns_megaplot.jpeg", width = 7,
#        height = 10, units = 'in')

# Compare site-level richness and covariate ----------------------
# Pull Zs from each item in list  
Zs <- lapply(biglist, function(x) x$BUGSoutput$sims.list$Z)

# Get avg. occurrence matrices  
Zs.mean <- lapply(Zs, apply, c(2,3), mean)

# Get site-level richness  
site.rich <- lapply(Zs.mean, rowSums)

# Convert list to data.frame
rich.frame <- as.data.frame(do.call(cbind, site.rich))
colnames(rich.frame) <- c('uninf', 'inf.weak', 'inf.mod', 
                          'inf.strong', 'misinf.weak', 
                          'misinf.mod', 'misinf.strong')

rich.frame$True <- colSums(tru)
rich.frame$Obs <- rowSums(apply(obs.data, c(1,3), max))
rich.frame$Cov <- cov

rich.long <- rich.frame %>%
  pivot_longer(uninf:Obs, names_to = 'model', 
               values_to = 'Richness') %>%
  mutate(mod.type = case_when(startsWith(model,"un")~"Uninformed",
                            startsWith(model,"mis")~"Misinformed",
                            startsWith(model,"inf")~"Informed",
                            model == "True"~"True",
                            model == "Obs"~"Observed"))

order <- c("Observed", "True", "Uninformed", "Informed", "Misinformed")
rich.long$mod.type <- factor(rich.long$mod.type, levels = order)
  
rich.plot <- ggplot(data = rich.long, aes(x = Cov, y = Richness,
                                         color = mod.type))+
  geom_point()+
  geom_smooth(aes(fill = mod.type), method = 'lm', alpha = 0.2)+
  labs(x = "Covariate")+
  expand_limits(y = 0)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(), 
        legend.title = element_blank())

# ggsave(rich.plot, filename = 'richplot.jpeg')

# Compare covariate responses ----------------------
get.cov <- function(jag){
  # Extract covariate estimates from jags object
  a1s <- jag$BUGSoutput$sims.list$a1
  
  a1s <- as.data.frame(a1s[,21:22])
  
  specnames <- logical()
  for(i in 21:22){
    specnames[i-20] <- paste("Spec", i, sep = "")
  }
  
  colnames(a1s) <- specnames

  # Pivot data frame for plotting
  a1.long <- a1s %>%
    pivot_longer(cols = everything(), names_to = "Spec", 
                 values_to = "a1")
  
  a1.long$Spec <- factor(a1.long$Spec, levels = specnames)

  a1.stat <- a1.long %>%
    group_by(Spec) %>%
    summarise(mean = mean(a1), lo = quantile(a1, 0.025), 
              hi = quantile(a1, 0.975)) %>%
    mutate(tru.resp = resp2cov[21:22])

  # Make interval plot
  plot <- ggplot(data = a1.stat, aes(x = Spec, y = mean))+
    geom_point(size = 1.5)+
    geom_errorbar(ymin = a1.stat$lo, ymax = a1.stat$hi, 
                  size = 1, width = 0.2)+
    geom_point(aes(y = tru.resp), color = "red", size = 1.5)+
    geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
    scale_y_continuous(limits = c(-10, 10))+
    labs(x = "Species", y = "Coefficient")+
    theme_bw(base_size = 14)+
    theme(panel.grid = element_blank())
  
  return(plot)
}

cov.plots <- lapply(biglist, get.cov)

plot.uninf <- cov.plots[[1]]
plot.inf <- cov.plots[[2]]/cov.plots[[3]]/cov.plots[[4]]
plot.misinf <- cov.plots[[5]]/cov.plots[[6]]/cov.plots[[7]]

plot.uninf/(plot.inf|plot.misinf)+
  plot_layout(heights = c(1,4))

# ggsave(allthecovs, filename = "allthecovs.jpeg", height = 10,
#        width = 8, units = "in")

# Compare typical bias of each prior method ---------------------
# Get series of site-level estimates from Zs
rich.bias <- function(jag){
  # Get site-level richness
  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.replicates <- apply(Zs, c(1,2), sum)

  # Subtract each replicate vector from the true value
  rich.diffs <- apply(Zs.replicates, 1, function(x) x-colSums(tru))
  
  # Convert to vector
  diffs.vec <- as.vector(rich.diffs)

  # Create histograms
  plot <- ggplot()+
      geom_bar(aes(x = diffs.vec), fill = 'darkgray')+
      geom_vline(xintercept = 0)+
      labs(x = "Richness")+
      theme_bw(base_size = 18)+
      theme(axis.title.y = element_blank(), 
            panel.grid = element_blank())
  
  return(plot)
}

rich.hists <- lapply(biglist, rich.bias)

plot.uninf <- plot_spacer()+rich.hists[[1]]+plot_spacer()+
  plot_layout(widths = c(1,2,1))
rich.hists.inf <- rich.hists[[2]]/rich.hists[[3]]/rich.hists[[4]]
rich.hists.misinf <- rich.hists[[5]]/rich.hists[[6]]/rich.hists[[7]]

allthehists <- plot.uninf/(rich.hists.inf|rich.hists.misinf)+
  plot_layout(heights = c(1,5))

# ggsave(allthehists, file = "allthehists.jpeg", height = 10,
#        width = 8, units = 'in')

# Function looking at undetected spec error --------------------
error.raster <- function(jag){
  Zs <- jag$BUGSoutput$sims.list$Z
  
  spec.names <- c("Spec21", "Spec22")

  Zs.mean <- data.frame(apply(Zs, c(2,3), mean)[,21:22])
  colnames(Zs.mean) <- spec.names
  Zs.mean$Site <- c(1:30)

  tru.frame <-as.data.frame(t(tru))[,21:22]
  colnames(tru.frame) <- spec.names
  tru.frame$Site <- c(1:30)

  tru.frame %>%
    gather("Spec21":"Spec22", key = "Species", value = "Occ") %>%
    {. ->> tru.frame}

  Zs.mean %>%
    gather("Spec21":"Spec22", key = "Species", value = "Occ") %>%
    left_join(tru.frame, by = c("Site", "Species")) %>%
    mutate(Error = Occ.y - Occ.x) %>%
    {. ->> merged.frame}

  merged.frame$Species = factor(merged.frame$Species, 
                                levels = spec.names)

  rast <- ggplot(data = merged.frame, 
                 aes(x = Species, y = Site, fill = Error))+
    geom_tile()+
    scale_fill_distiller(type = "div", 
                         limit = max(abs(merged.frame$Error)) *
                           c(-1, 1))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme_bw(base_size = 18)
  
  return(rast)
}

error.out <- lapply(biglist, error.raster)

plot.uninf <- plot_spacer() + error.out[[1]]+ plot_spacer()+
  plot_layout(widths = c(1,2,1))

plot.inf <- (error.out[[2]]/error.out[[3]]/error.out[[4]])&
  guides(fill = "none")

plot.misinf <- (error.out[[5]]/error.out[[6]]/error.out[[7]])&
  guides(fill = "none")

error.plots <- plot.uninf/(plot.inf|plot.misinf)+
  plot_layout(heights = c(1,5), guides = "collect")

# ggsave(error.plots, file = "error_plots.jpeg", height = 10, 
#        width = 8, units = "in")
