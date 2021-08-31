######################################################
# Appendix S3                                        #
# Code for data simulation and JAGS model            #
# Summer 2021                                        #
# E.M. Beasley                                       #
######################################################

# This model is written for R and uses the package R2jags
# to interface with the MCMC sampler JAGS

# The base model with uninformed priors is used to estimate species-
# specific occupancy and detection, and uses data augmentation to 
# estimate the number of hypothetical species that may have been undetected

# Models with aggregated priors introduce ecological information about
# these undetected species to the model, yielding more accurate 
# estimates and making species identity explicit

# Setup -----------------------------------------------
# Load packages
library(R2jags)
library(tidyverse)
library(MASS)
library(abind)
library(boot)

# Set seed
set.seed(39)

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

resp2cov <- sample(resp2cov) # reshuffle responses

# Covariate values for sites
cov <- sort(rnorm(n = nsite))

# Simulate occupancy data -------------------------------------
# Get species-level occupancy probs from a beta distribution
sim.occ <- rbeta(n = nspec+nmiss, shape1 = 2, shape2 = 4)
# Look at 95% interval
qbeta(p = c(0.025, 0.975), shape1 = 2, shape2 = 4)

# Write function to simulate true occupancy state
tru.mats <- function(spec=nspec+nmiss, site=nsite, alpha1=resp2cov){
  #Get site-level psi to account for covariates
  alpha0 <- logit(sim.occ)
  
  logit.psi <- matrix(NA, nrow = spec, ncol = site)
  
  for(i in 1:spec){
    logit.psi[i,] <- alpha0[i] + alpha1[i]*cov
  }
  
  psi <- inv.logit(logit.psi)
  
  #create list of abundance vectors
  nlist<-list()
  for(a in 1:spec){
    nlist[[a]] <- rbinom(n = site, size = 1, prob = psi[a,])
  }
  
  #Turn abundance vectors into abundance matrix
  ns<-do.call(rbind, nlist)
  
  return(ns)
}

tru <- tru.mats()

# Simulate detection process ---------------------------------
# Generate mean detection probabilities from beta dist
mean.p <- rbeta(n = nspec+nmiss, shape1 = 2, shape2 = 8)
# sorting by detection prob makes plotting easier
mean.p <- sort(mean.p, decreasing = T)
# Look at 95% interval
qbeta(p = c(0.025, 0.975), shape1=2, shape2=8)

# Generate detection histories
get.obs <- function(mat, specs){
  #Detection intercept and cov responses
  beta0<-logit(mean.p) #put it on logit scale

  #Logit link function
  logit.p <- array(NA, dim = c(nsite, nsurvey, specs))
  for(i in 1:specs){
    for(j in 1:nsite){
      for(k in 1:nsurvey){
        logit.p[j,,i] <- beta0[i] #detection covs would go here
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


  return(obsdata)
}

obs.data <- get.obs(mat = tru, specs = nspec+nmiss)

#Checkpoint: get observed data matrix
maxobs <- apply(obs.data, c(1,3), max)
colSums(maxobs)

# Remove species with no detections 
obs.data <- obs.data[,,-which(colSums(maxobs)==0)]

# Function to reorder values
# to put undetected species last
# this makes figures more clear
reorder <- function(x){
  if (length(dim(x)) == 0){
    nondets <- which(colSums(maxobs) == 0)
    copy <- x[nondets]
    x <- x[-nondets]
    new <- c(x, copy)
    return(new)
    }
  else {
    nondets <- which(colSums(maxobs) == 0)
    copy <- x[nondets,]
    x <- x[-nondets,]
    new <- rbind(x, copy)
    return(new)
    }
}

sim.occ <- reorder(sim.occ)
mean.p <- reorder(mean.p)

resp2cov <- reorder(resp2cov)

tru <- reorder(tru)

# Augment the observed dataset ---------------------------------
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

            inf.mean0 <- c(0, -2,0, 0)
            inf.mean1 <- c(0, 0,-3, 0)
            
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

            inf.mean0 <- c(0, -2,0, 0)
            inf.mean1 <- c(0, 0,-3, 0)
            
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

            inf.mean0 <- c(0, -2,0, 0)
            inf.mean1 <- c(0, 0,-3, 0)
            
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

            inf.mean0 <- c(0, 2,0, 0)
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

            inf.mean0 <- c(0, 2,0, 0)
            inf.mean1 <- c(0, 0,3, 0)
            
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

            inf.mean0 <- c(0, 2,0, 0)
            inf.mean1 <- c(0, 0,3, 0)
            
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

# Write Model ----------------------------
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

# Write function for sending model to gibbs sampler --------------
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

# these are commented out because running them can take a while

# mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug,
#                         cov = cov, spec = nspec,
#                         textdoc = 'aug_model.txt', aug = nmiss,
#                         burn = 2500, iter = 10000, thin = 10)
# saveRDS(mod.uninf, file = "mod_uninf.rds")
# 
# inf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug,
#                        cov = cov, spec = nspec,
#                        textdoc = 'aug_model.txt', aug = nmiss,
#                        priors = weakinf, burn = 3000,
#                        iter = 10000, thin = 5)
# saveRDS(inf.weak, file = "inf_weak.rds")
# 
# inf.mod <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug,
#                       cov = cov, spec = nspec,
#                       textdoc = 'aug_model.txt', aug = nmiss,
#                       priors = modinf, burn = 8000,
#                       iter = 12000, thin = 3)
# saveRDS(inf.mod, file = "inf_mod.rds")
# 
# inf.strong <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug,
#                          cov = cov, spec = nspec,
#                          textdoc = 'aug_model.txt', aug = nmiss,
#                          priors = stronginf, burn = 8000,
#                          iter = 12000, thin = 3)
# saveRDS(inf.strong, file = "inf_strong.rds")
# 
# misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, 
#                           cov = cov, spec = nspec, 
#                           textdoc = 'aug_model.txt',
#                           aug = nmiss, priors = weakmisinf, 
#                           burn = 5000, iter = 10000, thin = 5)
# saveRDS(misinf.weak, file = "misinf_weak.rds")
# 
# misinf.mod <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, 
#                          cov = cov, spec = nspec, 
#                          textdoc = 'aug_model.txt', aug = nmiss,
#                          priors = modmisinf, burn = 2500,
#                          iter = 10000, thin = 10)
# saveRDS(misinf.mod, file = "misinf_mod.rds")
# 
# misinf.strong <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, 
#                             cov = cov, spec = nspec, 
#                             textdoc = 'aug_model.txt',
#                             aug = nmiss, priors = strongmisinf,
#                             burn = 2000, iter = 10000, thin = 5)
# saveRDS(misinf.strong, file = "misinf_strong.rds")

# Load models -------------------------------
uninf <- readRDS("mod_uninf.rds")
inf.weak <- readRDS("inf_weak.rds")
inf.mod <- readRDS("inf_mod.rds")
inf.strong <- readRDS("inf_strong.rds")
misinf.weak <- readRDS("misinf_weak.rds")
misinf.mod <- readRDS("misinf_mod.rds")
misinf.strong <- readRDS("misinf_strong.rds")

