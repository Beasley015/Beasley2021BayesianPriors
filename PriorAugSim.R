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
library(gridExtra)
library(grid)
library(agricolae)
library(ggridges)
library(data.table)

# Set seed
set.seed(15)

# Global variables
nspec <- 20
nmiss <- 2 # Species present but not detected during sampling
naug <- 3 # Species never detected; used to set prior on N
nsite <- 30
nsurvey <- 4

Ks <- rep(nsurvey, nsite) # vector of surveys per site

# Write JAGS scripts for priors on undetected species
# Uninformative priors
uninf <- "for(i in 1:(spec+aug)){
          #Create priors from hyperpriors
            w[i] ~ dbern(omega)
            
            a0[i] ~ dnorm(a0.mean, tau.a0)
            a1[i] ~ dnorm(a1.mean, tau.a1)

            b0[i] ~ dnorm(b0.mean, tau.b0)"


# Weakly informative priors
weakinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
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

# Moderately informative
modinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
                            
            inf.mean1 <- c(0, 0, -3, 0)
            
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

# Strongly informative
stronginf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 0, -3, 0)
            
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

# Weakly mis-specified
weakmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, 3, 0)
            
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

# Moderately mis-specified
modmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, 3, 0)
            
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

# Strongly mis-specified
strongmisinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, 3, 0)
            
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

# Function for simulated covs and coefficients -------------------
sim.covs <- function(){
  # Vector of covariate responses (detected species)
  resp2cov <- c(rnorm(n = 6, sd = 0.25),
                rnorm(n = 7, mean = 3, sd = 0.25),
                rnorm(n = 7, mean = -3, sd = 0.25))

  # Add undetected species
  resp2cov[21:22] <- c(rnorm(n = 1, mean = 0, sd = 0.25),
                       rnorm(n = 1, mean = -3, sd = 0.25))

  # Covariate values for sites
  cov <- sort(rnorm(n = nsite))
  
  covs.out <- list(resp2cov, cov)
  
  return(covs.out)
}

# Function to simulate occupancy data ---------------------------
occ.func <- function(resp2cov, cov){
  # Get occupancy probs from a beta distribution
  sim.occ <- rbeta(n = nspec, shape1 = 2, shape2 = 4)
  # Keep undetected species consistent in all sims
  sim.occ[21:22] <- c(0.1, 0.4)

  #Get site-level occupancy prob (psi) to account for covariates
  alpha0 <- logit(sim.occ)
  
  logit.psi <- matrix(NA, nrow = nspec+nmiss, ncol = nsite)
  
  for(i in 1:(nspec+nmiss)){
    logit.psi[i,] <- alpha0[i] + resp2cov[i]*cov
  }
  
  psi <- inv.logit(logit.psi)
  
  #create list of occupancy vectors for each species
  nlist<-list()
  for(a in 1:(nspec+nmiss)){
    nlist[[a]] <- rbinom(n = nsite, size = 1, prob = psi[a,])
  }
  
  #Turn vectors into presence/absence matrix
  ns<-do.call(rbind, nlist)
  
  # Fill any 0s so all species occupy at least 1 site
  if(any(rowSums(ns)==0)){
    missing <- which(rowSums(ns)==0)
    
    for(i in 1:length(missing)){
      ns[missing[i],
         psi[missing[i],
             which(psi[missing[i],]==max(psi[missing[i],]))]]<- 1
    }
  }
  return(list(sim.occ, ns, psi))
}

# Function to simulate detection process -------------------------
det.func <- function(mat, psi){
  # Generate mean detection probabilities from beta dist
  mean.p <- rbeta(n = nspec, shape1 = 2, shape2 = 8)
  mean.p <- sort(mean.p, decreasing = T) # sorting keeps it tidy

  # add undetected species
  mean.p[21:22] <- 0

  #Detection intercept and cov responses
  beta0<-logit(mean.p) #put it on logit scale

  #Logit link function: accounts for site, survey covariates
  # if covariates are present
  logit.p <- array(NA, dim = c(nsite, nsurvey, (nspec+nmiss)))
  for(i in 1:(nspec+nmiss)){
    for(j in 1:nsite){
      for(k in 1:nsurvey){
        logit.p[j,,i] <- beta0[i] # add detection covariates here
      }
    }
  }

  p <- plogis(logit.p)

  #Simulate observation data
  L<-list()

  for(b in 1:(nspec+nmiss)){
    y<-matrix(NA, ncol = nsite, nrow = nsurvey)
    for(a in 1:nsurvey){
      y[a,]<-rbinom(n = nsite, size = 1, prob = p[,,b]*mat[b,])
    }
    L[[b]]<-t(y)
  }

  #Smash it into array
  obsdata <- array(as.numeric(unlist(L)), dim=c(nsite, nsurvey, 
                                              (nspec+nmiss)))

  #Get site-level p
  ps <- apply(p, c(1,3), mean)

  # Get observed data matrix
  maxobs <- apply(obsdata, c(1,3), max)

  # Make sure first 20 species are detected
  if(any(colSums(maxobs[,1:20])==0)){
    missing <- which(colSums(maxobs[,1:20])==0)
  
    for(i in 1:length(missing)){
      obsdata[which(psi[missing[i],]==max(psi[missing[i],])),
               sample(1:nsurvey, size = 1),
               missing[i]] <- 1
    }
  
    maxobs <- apply(obsdata, c(1,3), max)
  }
  
  # Add augmented species - not present, puts prior on N
  augdata <- array(0, dim = c(nsite, nsurvey, naug))
  
  obsdata <- abind(obsdata, augdata, along = 3)
  
  return(obsdata) 
}

# Function to create all community data ------------------
comm.sim <- function(){
  # Run all simulation functions
  covs <- sim.covs()
  occ.sim <- occ.func(resp2cov = covs[[1]], cov = covs[[2]])
  det.sim <- det.func(mat = occ.sim[[2]], psi = occ.sim[[3]])
  
  # Save relevant parameters to list
  out.list <- list(cov = covs[[2]], sim.occ = occ.sim[[1]], 
                   obs = det.sim, tru = occ.sim[[2]])
  return(out.list)
}

# Function for sending model to gibbs sampler --------------
VivaLaMSOM <- function(J = nsite, K = nsurvey, obs, spec, aug = 0, 
                       cov, sim.occ, textdoc, 
                       priors = uninf, burn = 2500, iter = 8000, 
                       thin = 10){
  
  # Write model for augmented datasets
  if(textdoc == 'aug_model.txt')
    write.model(priors = priors)
  
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, cov = cov,
                   sim.occ = sim.occ)
  if(textdoc == 'aug_model.txt'){
    datalist$aug <- aug
  }
  
  # Specify parameters
  parms <- c('a0','a1','Z','N','a1.mean')
  
  # Initial values
  maxobs <- apply(obs, c(1,3), max) # starting w observed data
  # helps avoid parent node incompatibility
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
  model <- jags(model.file = textdoc, data = datalist, 
                n.chains = 3, parameters.to.save = parms, 
                inits = init.values, n.burnin = burn, 
                n.iter = iter, n.thin = thin)
  
  return(model)
}

# Write JAGS scripts ----------------------------
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
    
        #Estimate detection of i at point j during day k
        for(k in 1:K[j]){
          logit(p[j,k,i]) <-  b0[i]
          mu.p[j,k,i] <- p[j,k,i]*Z[j,i] 
          #The addition of Z means that detecting a species depends on its occupancy
          obs[j,k,i] ~ dbern(mu.p[j,k,i])
    }
    }
    }
    
    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
    n0<-w[(spec+1):(spec+aug)]
    N<-spec+sum(n0)
    }
    ")
 writeLines(mod, "aug_model.txt") 
}

# Function to fit models ------------------------------------
fit.mods <- function(cov, obs, sim.occ){
  # Non-Augmented model just in case
  # mod.noaug <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
  #                         textdoc = 'noaug.txt')

  # Mod with uninformative priors
  mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                          textdoc = 'aug_model.txt', cov = cov,
                          obs = obs, sim.occ = sim.occ,
                          aug = nmiss+naug, burn = 2500, 
                          iter = 10000, thin = 10)

  # Mods with informative priors of varying strength
  inf.weak <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                         textdoc = 'aug_model.txt', cov = cov,
                         obs = obs, sim.occ = sim.occ,
                         aug = nmiss+naug, priors = weakinf, 
                         burn = 3000, iter = 10000, thin = 5)

  inf.mod <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                        textdoc = 'aug_model.txt', cov = cov,
                        obs = obs, sim.occ = sim.occ,
                        aug = nmiss+naug, priors = modinf, 
                        burn = 8000, iter = 12000, thin = 3)

  inf.strong <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                           textdoc = 'aug_model.txt', cov = cov,
                           obs = obs, sim.occ = sim.occ, 
                           aug = nmiss+naug, priors = stronginf,
                           burn = 8000, iter = 12000, thin = 3)

  # Mods with mis-specified priors of varying strength
  misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                            textdoc = 'aug_model.txt', cov = cov,
                            obs = obs, sim.occ = sim.occ,
                            aug = nmiss+naug, priors = weakmisinf,
                            burn = 5000, iter = 10000, thin = 5)

  misinf.mod <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                           textdoc = 'aug_model.txt', cov = cov,
                           obs = obs, sim.occ = sim.occ,
                           aug = nmiss+naug, priors = modmisinf,
                           burn = 2500, iter = 10000, thin = 10)

  misinf.strong <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                              textdoc = 'aug_model.txt', cov = cov,
                              obs = obs, sim.occ = sim.occ, 
                              aug = nmiss+naug, 
                              priors = strongmisinf, burn = 2000, 
                              iter = 10000, thin = 5)
  
  out.list <- list(mod.uninf, inf.weak, inf.mod, inf.strong,
                   misinf.weak, misinf.mod, misinf.strong)
  
  return(out.list)
}

# Run that sucker ----------------------
forth.eorlingas <- function(iters){
  # List for community sim outputs
  sim.res <- list()
  
  # Run the models
  for(i in 1:iters){
    sim.outs <- comm.sim()
    sim.res <- append(sim.res, sim.outs)
    
    results <- fit.mods(cov = sim.outs[[1]], obs = sim.outs[[3]],
                        sim.occ = sim.outs[[2]])
    
    filename <- paste("./Outputs/out", i, ".rds", sep = "")
    
    saveRDS(results, file = filename)
    
    print(i)
  }
  
  # save sim outputs
  saveRDS(sim.res, file = "simres.rds")
}

# forth.eorlingas(iters = 50) # warning: will take several hours
# (approx. 12-15) unless you set up parallel processing

# there is also a warning I haven't fixed 
# It will not affect final results

# Function to load results ---------------
get.outs <- function(param){
  # Pull file names from outputs folder
  filenames <- list.files("./Outputs")
  
  # Create blank list to store results
  param.list <- list()
  # Vector of model types
  modtype <- c('mod.uninf', 'inf.weak', 'inf.mod', 'inf.strong',
               'misinf.weak', 'misinf.mod', 'misinf.strong')
  
  # Pulls desired results from each model
  for(i in 1:length(filenames)){
    out <- readRDS(file = paste("./Outputs/", filenames[i],
                                sep = ""))
    
    for(j in 1:length(modtype)){
      param.out <- out[[j]]$BUGSoutput$sims.list[param]
      
      if(length(param.list) < length(modtype)){
        param.list[[j]] <- param.out
      } else{
        param.list[[j]] <- append(param.list[[j]], param.out)
      }
    }
    
    print(paste("i = ", i, sep = ""))
  }
  
  # add model names to results output
  names(param.list) <- modtype
  
  return(param.list)
}

# Compare hyperparameter distribution --------------------
hyper.means <- get.outs(param = "a1.mean")

fewer.means <- hyper.means[c("mod.uninf", "inf.strong", 
                             "misinf.strong")]
names(fewer.means) <- c("Uninformative", "StrongInformative",
                        "StrongMisSpecified")

for(i in 1:length(fewer.means)){
  fewer.means[[i]] <- lapply(fewer.means[[i]],
                             function(x) 
                               data.frame(cbind(as.vector(x), 
                                     names(fewer.means)[[i]])))
}

fewer.means <- do.call('c', fewer.means)

reps <- rep(1:50, 3)
for(i in 1:length(fewer.means)){
  fewer.means[[i]] <- cbind(fewer.means[[i]], reps[i])
}

means.frame <- do.call(rbind, fewer.means)
means.frame$X1 <- as.numeric(means.frame$X1)
colnames(means.frame) <- c("a1.mean", "modtype", "rep")

ggplot(data = means.frame, aes(x = a1.mean, 
                               group = interaction(modtype,rep), 
                               color = modtype))+
  geom_density()+
  xlim(c(-4, 4))+
  labs(x = "a1 Mean", y = "Density")+
  scale_color_viridis_d(name = "Model Type", end = 0.9)+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename="a1posteriors.jpeg", width = 4, height = 4, 
#        units = "in")

ggplot(data = means.frame[means.frame$rep %in% 1:9,], 
       aes(x = a1.mean, color = modtype))+
  geom_density()+
  facet_wrap(vars(as.character(rep)))+
  labs(x = "Community Mean (a1)", y = "Density")+
  scale_color_viridis_d(name = "Model Type", end = 0.9)+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "a1posteriors_facet.jpeg", width = 6, height = 6,
#        units = "in")

# Compare Ns -------------------------
# Load all n's into R
all.n <- get.outs(param = "N")

# Calculate measure of centrality
ns.center <- function(jag, center = "mode"){
  if(center == "mode"){
    # Create list for storage of tables
    jagtab <- list()
    # For loop to calculate mode for each element
    for(i in 1:50){
      #Convert to character
      jag[[i]] <- as.character(jag[[i]])
      # Create table
      jagtab[[i]] <- table(jag[[i]])
      
      # Find most common value in table
      jag[[i]] <- names(table(jag[[i]]))[which.max(jagtab[[i]])]
    }
    # Change to vector of modes
    x = unlist(jag)
    return(x)
    
  } else if(center == "median"){
    for(i in 1:50){
    # Get median
     jag[[i]] <- median(as.vector(jag[[i]]))
    }
    # Change to vector
    x = unlist(jag)
    return(x)
    
  } else if(center == "mean"){
    for(i in 1:50){
    # Get mean
    jag[[i]] <- mean(jag[[i]])
    }
    # Change to vector
    x = unlist(jag)
    return(x)
  }
}

# Get summary of metacommunity richness using different
# measures of centrality 
ns.means <- lapply(all.n, ns.center, center = "mean")
ns.medians <- lapply(all.n, ns.center, center = "median")
ns.modes <- lapply(all.n, ns.center, center = "mode")

# Plot it
ns.plot <- function(dat, center){
  if(center == "mean"){
    Ns.plot <- ggplot(data = data.frame(ns = dat), aes(x = ns))+
      geom_histogram(binwidth = 0.5, color = 'lightgray')+
      geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"),
                 size = 1.5)+
      scale_y_continuous(expand = c(0,0))+
      theme_classic(base_size = 12)+
      theme(axis.text.y = element_blank(), 
            axis.title = element_blank(), 
            plot.margin = unit(c(0,0,0,0), units = "point"),
            legend.position = "None")
  } else{
  ns.frame <- dat %>%
    table() %>%
    data.frame() %>%
    {. ->> ns.frame}
  colnames(ns.frame) <- c("N_Species", "Freq")

  
  Ns.plot <- ggplot(data = ns.frame, 
                    aes(x = as.integer(as.character(N_Species)), 
                        y = Freq))+
    geom_col(color = 'lightgray')+
    geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"),
               size = 1.5)+
    scale_y_continuous(expand = c(0,0))+
    theme_classic(base_size = 12)+
    theme(axis.text.y = element_blank(), 
          axis.title = element_blank(), 
          plot.margin = unit(c(0,0,0,0), units = "point"),
          legend.position = "None")
  
  }
  return(Ns.plot)
}

# Look at summarised results
nouts.mode <- lapply(ns.modes, ns.plot, center = "mode")

# Define plot layout
layout <- "
#AA#
BBEE
CCFF
DDGG
"

histos <- nouts.mode

# Create figure
Ns.base <- histos[[1]]+histos[[2]]+histos[[3]]+histos[[4]]+
  histos[[5]]+histos[[6]]+histos[[7]]+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "a")+
  plot_layout(guides = "collect")&
  xlim(19.5, 25.5)

# Add labels
gn <- patchworkGrob(Ns.base)
Ns.megaplot <- grid.arrange(gn, bottom = textGrob("Species Richness (N)", 
                                                  gp=gpar(fontsize = 14), hjust = 0.8))

# Save figure
# ggsave(Ns.megaplot, filename = "ns_mode.jpeg", width = 6,
#        height = 5, units = 'in', dpi = 600)

# Alternative: Ridge plot
# Assign reps to each
list.combined <- lapply(all.n, function(y) 
  as.data.frame(do.call("rbind", y)))

list.combined2 <- lapply(list.combined, function(y)
  data.frame(y, rep = rep(1:50, each = nrow(y)/50)))

for(i in 1:length(list.combined2)){
  list.combined2[[i]]$model <- names(list.combined2)[i]
}

ns.frame <- rbindlist(list.combined2)

ns.frame <- ns.frame %>%
  group_by(V1, rep, model) %>%
  summarise(count = n())

test <- filter(ns.frame, model == "misinf.strong")

strongmis.ridge <- ggplot(data = test, aes(x = V1, y = rep, 
                                         group = rep))+
  geom_ridgeline(aes(height = count, scale = 0.1), alpha = 0.5)+
  geom_vline(xintercept = 22, linetype = "dashed")+
  lims(y = c(0,350))+
  labs(x = "Regional Richness (Estimated)", y = "Frequency")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave(strongmis.ridge, filename = "strongmis_ridge.jpeg", width = 3, 
       height = 3, units = "in")

(uninf.ridge + weakmis.ridge)/
  (modmis.ridge + strongmis.ridge)&
  plot_annotation(tag_levels = "a")

# ggsave(filename = "misplotsN.jpeg", width = 6, height = 5, 
#        units = "in")

# Compare covariate responses ----------------------
# Load in cov responses
all.a1 <- get.outs(param = "a1")

# Reformat lists
list.to.frame <- function(x){
  # Combine all iterations of each model
  list.combined <- lapply(x, function(y) 
    as.data.frame(do.call("rbind", y)))
  
  # Create species list
  specs <- logical()
  for(i in 1:(nspec+naug+nmiss)){
    specs[i] <- paste("Spec", i, sep = "")
  }
  
  # rename columns
  list.cols <- lapply(list.combined, setNames, specs)
  
  # Add column to each data frame
  modnames <- names(list.cols)
  list.named <- mapply(cbind, list.cols, "model" = modnames, 
                     SIMPLIFY = F)
  
  # Merge it all into data frame
  big.ass.data.frame <- do.call(rbind, list.named)
  
  return(big.ass.data.frame)
}

chonk <- list.to.frame(all.a1)

# pivot to long format
a1s <- chonk %>%
    pivot_longer(cols = -model, names_to = "Species",
                 values_to = "a1")
  
# Calculate means and quantiles
a1.stat <- a1s %>%
  group_by(Species, model) %>%
  summarise(mean = mean(a1), lo = quantile(a1, 0.025), 
              hi = quantile(a1, 0.975)) %>%
  mutate(Type = case_when(startsWith(as.character(model), "mod") ~ 
                            "Uninformative", 
                          startsWith(as.character(model), "inf") ~
                            "Informative", 
                          startsWith(as.character(model), "mis") ~
                            "Mis-specified")) %>%
  mutate(Type = factor(Type, levels = c("Uninformative", 
                                        "Informative",
                                        "Mis-specified"))) %>%
  mutate(Weight = case_when(endsWith(as.character(model), "weak") ~
                              "Weak",
                            endsWith(as.character(model), "mod") ~
                              "Moderate",
                            endsWith(as.character(model), "strong") ~
                              "Strong")) %>%
  mutate(Weight = factor(Weight, levels = c("Weak", "Moderate", 
                                            "Strong")))
  
smol.a1.21 <- filter(a1.stat, Species == "Spec21")
smol.a1.22 <- filter(a1.stat, Species == "Spec22")

# Make interval plots
covplot.21 <- ggplot(data = smol.a1.21, aes(x = Weight, y = mean))+
  geom_point(size = 1.5)+
  geom_errorbar(ymin = smol.a1.21$lo, ymax = smol.a1.21$hi,
                  size = 1, width = 0.2)+
  geom_point(aes(y = 0), color = "red", size = 1.5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0))+
  labs(y = "Coefficient")+
  theme_bw(base_size = 14)+
  facet_grid(~Type, scales = "free_x", space = "free_x",
             switch = "x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(), axis.title.x = element_blank())

covplot.22 <- ggplot(data = smol.a1.22, aes(x = Weight, y = mean))+
  geom_point(size = 1.5)+
  geom_errorbar(ymin = smol.a1.22$lo, ymax = smol.a1.22$hi,
                size = 1, width = 0.2)+
  geom_point(aes(y = -3), color = "red", size = 1.5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0))+
  labs(x = "Model", y = "Coefficient")+
  theme_bw(base_size = 14)+
  facet_grid(~Type, scales = "free_x", space = "free_x",
             switch = "x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(), axis.title.x = element_blank())

covs <- (covplot.21/covplot.22)+
  plot_annotation(tag_levels = "a")

# ggsave(covs, filename = "undet_cov.jpeg", dpi = 600, width = 8,
#        height = 6, units = "in")

# Compare site-level richness btn models ----------------------
# Get true values
sim.res <- read_rds(file = "simres.rds")
indices <- seq(from = 4, to = 200, by = 4)

tru <- list()
for(i in 1:length(indices)){
  tru[[i]] <- sim.res[[indices[i]]]
}

# Get differences between true and est values
get.diff <- function(){
  # Load in estimated vals
  zs <- get.outs(param = "Z")
  
  z.avg <- list()
  for(i in 1:7){
    z.avg[[i]] <- list()
    for(j in 1:50){
      z.avg[[i]][[j]] <- apply(zs[[i]][[j]], c(2,3), mean)[,-c(23:24)]
    }
  }

  # Get species richness
  tru.rich <- lapply(tru, colSums)

  z.rich <- list()
  for(i in 1:length(z.avg)){
      z.rich[[i]] <- lapply(z.avg[[i]], rowSums)
  }

  # subtract true and estimated richness
  z.site <- list()
  for(i in 1:length(z.rich)){
    z.site[[i]] <- lapply(z.rich[[i]], 
                          function(x) cbind(x, tru.rich[[i]]))
  }
  names(z.site) <- names(zs)
  
  return(z.site)
}

z.site <- get.diff()

# Get true coefficients
indices <- seq(from = 1, to = length(sim.res), by = 4)
coefs <- list()
for(i in 1: length(sim.res)){
  coefs[[i]] <- sim.res[[indices[i]]]
}

# Create vectors of site, rep names
sitenames <- logical()
for(i in 1:nsite){
  sitenames[i] <- paste("Site", i, sep = "")
}

repnames <- logical()
for(i in 1:50){
  repnames[i] <- paste("Rep", i, sep = "")
}

modnames <- names(z.site)

# Coerce all to data frames
coef.frame <- as.data.frame(do.call(cbind, coefs))
colnames(coef.frame) <- repnames
coef.frame$Site <- sitenames
coef.frame <- pivot_longer(coef.frame, -Site, names_to = "Rep", 
                           values_to = "Cov")

# More shenanigans to get the data right
diff.frame <- lapply(z.site, function(x) as.data.frame(do.call(rbind, x)))
diff.frame <- do.call(rbind, diff.frame)
diff.frame$Model <- rep(modnames, each = nrow(diff.frame)/7)
colnames(diff.frame) <- c("Est", "Tru", "Model")
diff.frame$Site <- rep(sitenames, 7)
diff.frame$Rep <- rep(repnames, each = 30)

# the final data frame
big.ass.frame <- left_join(diff.frame, coef.frame, 
                           by = c("Rep", "Site")) %>%
  mutate(Model = factor(Model, levels = unique(diff.frame$Model))) %>%
  mutate(Diff = (Tru-Est)/Tru) %>%
  group_by(Model, Rep) %>%
  summarise(mean.diff = mean(Diff), med.diff = median(Diff))
# using median because skewed low

# Use aov to get effect sizes, even though aov isn't used
med.mod <- aov(data = big.ass.frame, med.diff~Model)
mean.mod <- aov(data = big.ass.frame, mean.diff~Model)

# Effect sizes
summary(mean.mod)[[1]]$`Sum Sq`[1]/sum(summary(mean.mod)[[1]]$`Sum Sq`)

# use HSD to assign groups
med.hsd <- HSD.test(y = med.mod, trt = "Model")
mean.hsd <- HSD.test(y = mean.mod, trt = "Model")

hsd.df <- data.frame(Model = rownames(mean.hsd$groups),
                     group = mean.hsd$groups$groups)

# More manipulation to get labels right on figure
big.ass.frame <- left_join(big.ass.frame, hsd.df, by = "Model") %>%
  mutate(Model = factor(Model, 
                        levels = unique(big.ass.frame$Model))) %>%
  mutate(Type = case_when(startsWith(as.character(Model), "mod") ~
                            "Uninformative",
                          startsWith(as.character(Model), "inf") ~
                            "Informative",
                          startsWith(as.character(Model), "mis") ~
                            "Mis-specified")) %>%
  mutate(Type = factor(Type, levels = unique(Type))) %>%
  mutate(Weight = case_when(endsWith(as.character(Model), "weak") ~
                              "Weak",
                            endsWith(as.character(Model), "mod") ~
                              "Moderate",
                            endsWith(as.character(Model), "strong") ~
                              "Strong")) %>%
  mutate(Weight = factor(Weight, levels = unique(Weight)))
  
# Going with box plots for now
ggplot(data = big.ass.frame, aes(x = Type, y = mean.diff,
                                 fill = Weight))+
  geom_boxplot(outlier.shape = NA, varwidth = T)+
  scale_fill_viridis_d(na.value = "lightgray")+
  geom_text(aes(label = group, y = 0.3), 
            position = position_dodge(width = 1))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  labs(y = "Mean Difference (Scaled)")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank())

# ggsave("siterich.jpeg", height = 4, width = 6, units = "in")

# Look at occ estimates for missing species ---------------------
# Get true values
sim.res <- list()
for(i in 1:50){
  sim.res[[i]] <- comm.sim()
}

tru <- list()
for(i in 1:length(sim.res)){
  tru[[i]] <- sim.res[[i]][[4]][c(21:22),]
}

# Get differences between true and est values
get.undet <- function(){
  # Load in estimated vals
  zs <- get.outs(param = "Z")
  
  z.avg <- list()
  for(i in 1:7){
    z.avg[[i]] <- list()
    for(j in 1:50){
      z.avg[[i]][[j]] <- apply(zs[[i]][[j]], c(2,3), mean)[,21:22]
    }
  }
  
  # subtract true and estimated richness
  z.diff <- list()
  for(i in 1:length(z.avg)){
    z.diff[[i]] <- lapply(z.avg[[i]], function(x) tru[[i]]-t(x))
  }
  names(z.diff) <- names(zs)
  
  return(z.diff)
}

undet.err <- get.undet()

# Add cols denoting model and species
df.list <- lapply(undet.err, function(x) do.call(rbind, x))
for(i in 1:length(df.list)){
  df.list[[i]] <- as.data.frame(df.list[[i]])
  df.list[[i]]$model <- names(df.list)[i]
  df.list[[i]]$species <- rep(c("Spec21", "Spec22"), 50)
}

# make it a data frame
undet.frame <- do.call(rbind, df.list)

# Create vectors of site names
sitenames <- logical()
for(i in 1:nsite){
  sitenames[i] <- paste("Site", i, sep = "")
}

colnames(undet.frame) <- c(sitenames, "model", "Species")

# Pull simulated covariates
cov.list <- lapply(sim.res, '[[', 1)
cov.vec <- do.call(c, cov.list)

# Pivot longer and add column for covariate
undet.long <- undet.frame %>%
  pivot_longer(cols = Site1:Site30, names_to = "Site", 
               values_to = "Occ") %>%
  mutate(Site = factor(.$Site, levels = sitenames)) %>%
  mutate(cov = rep(cov.vec, 14))

# Get data frame for each species
undet.s21 <- filter(undet.long, Species == "Spec21")
undet.s22 <- filter(undet.long, Species == "Spec22")

unique.mods <- unique(undet.s21$model)

undet.figs <- function(df, mods){
  # Create list to hold figures
  plt <- list()

  for(i in 1:length(mods)){
    # create data frame 
    small.df <- df[df$model == mods[i], ]
    
    # plot results of one model type
    plt[[i]] <- ggplot(data = small.df, aes(x = cov, y = Occ))+
      geom_point()+
      geom_hline(yintercept = 0, linetype = "dashed")+
      scale_y_continuous(limits = c(-1,1))+
      labs(x = "Covariate", y = "True Occupancy-Estimated Occupancy")+
      theme_bw()+
      theme(panel.grid = element_blank(), 
            axis.title = element_blank())
  }
  
  return(plt)
}

# Plots
s21 <- undet.figs(df = undet.s21, mods = unique.mods)
s22 <- undet.figs(df = undet.s22, mods = unique.mods)

layout <- "
#AA#
BBEE
CCFF
DDGG
"

# Put them all together
s21.full <- s21[[1]]+s21[[2]]+s21[[3]]+s21[[4]]+s21[[5]]+
  s21[[6]]+s21[[7]]+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "a")

gn <- patchworkGrob(s21.full)
s21final <- grid.arrange(gn, bottom = textGrob("Covariate", 
                              gp=gpar(fontsize = 16), hjust = 0.4),
                            left = textGrob("True-Estimated Occupancy"
                                          ,rot = 90,
                                          gp = gpar(fontsize = 16)))

# ggsave(s21final, filename = "s21occ.jpeg", width = 10,
#        height = 10, units = 'in', dpi = 600)

s22.full <- s22[[1]]+s22[[2]]+s22[[3]]+s22[[4]]+s22[[5]]+
  s22[[6]]+s22[[7]]+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "a")

gn <- patchworkGrob(s22.full)
s22final <- grid.arrange(gn, bottom = textGrob("Covariate", 
                                               gp=gpar(fontsize = 16),
                                               hjust = 0.4),
                         left = textGrob("True-Estimated Occupancy"
                                         ,rot = 90,
                                         gp = gpar(fontsize = 16)))

# ggsave(s22final, filename = "s22occ.jpeg", width = 10,
#        height = 10, units = 'in', dpi = 600)

# Model coefficients for detected/other augmented species ---------
all.a1 <- get.outs(param = "a1")

# Reformat lists
list.to.frame <- function(x){
  # Assign reps to each
  list.combined <- lapply(x, function(y) 
    as.data.frame(do.call("rbind", y)))
  
  list.combined2 <- lapply(list.combined, function(y)
    data.frame(y, rep = rep(1:50, each = nrow(y)/50)))
  
  # Create species list
  specs <- logical()
  for(i in 1:(nspec+naug+nmiss)){
    specs[i] <- paste("Spec", i, sep = "")
  }
  
  # rename columns
  list.cols <- lapply(list.combined2, setNames, c(specs, "rep"))
  
  # Add column to each data frame
  modnames <- names(list.cols)
  list.named <- mapply(cbind, list.cols, "model" = modnames, 
                       SIMPLIFY = F)
  
  # Merge it all into data frame
  big.ass.data.frame <- do.call(rbind, list.named)
  
  return(big.ass.data.frame)
}

chonk <- list.to.frame(all.a1)

# pivot to long format
a1s <- chonk %>%
  pivot_longer(cols = -c(model, rep), names_to = "Species",
               values_to = "a1")

# Calculate means and quantiles
a1.stat <- a1s %>%
  group_by(Species, model,rep) %>%
  summarise(mean = mean(a1), lo = quantile(a1, 0.025), 
            hi = quantile(a1, 0.975)) %>%
  mutate(Type = case_when(startsWith(as.character(model), "mod") ~ 
                            "Uninformative", 
                          startsWith(as.character(model), "inf") ~
                            "Informative", 
                          startsWith(as.character(model), "mis") ~
                            "Mis-specified")) %>%
  mutate(Type = factor(Type, levels = c("Uninformative", 
                                        "Informative",
                                        "Mis-specified"))) %>%
  mutate(Weight = case_when(endsWith(as.character(model), "weak") ~
                              "Weak",
                            endsWith(as.character(model), "mod") ~
                              "Moderate",
                            endsWith(as.character(model), "strong") ~
                              "Strong")) %>%
  mutate(Weight = factor(Weight, levels = c("Weak", "Moderate", 
                                            "Strong")))

spec.subset <- c(paste("Spec", sample(1:20, size = 3), sep = ""),
                 "Spec24", "Spec25")

smol.frame <- a1.stat %>%
  filter(Species %in% spec.subset) %>%
  filter(rep < 4) %>%
  filter((Weight != "Weak") %>% replace_na(TRUE)) %>%
  filter((Weight != "Moderate") %>% replace_na(TRUE)) %>%
  mutate(Species = factor(Species, levels = spec.subset))

noprior1 <- ggplot(data = smol.frame[smol.frame$rep == 1,], 
       aes(x = Species, y = mean,color = Type))+
  geom_point(position = position_dodge(width=.9))+
  geom_errorbar(aes(ymin = lo, ymax = hi), 
                position = position_dodge(width = .9))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_color_viridis_d(end = 0.9, name = "Model Type")+
  labs(y = "Coefficient Estimate (a1)")+
  theme_bw()+
  theme(panel.grid = element_blank())

# ggsave(filename = "NoPrior_rep3.jpeg", width = 4, height = 3,
#        units = "in")

(noprior1/noprior2/noprior3)+
  plot_layout(guides = "collect")&
  plot_annotation(tag_levels = "a")

# ggsave(filename = "NoPrior_allrep.jpeg", width = 5, height = 8,
#        units = "in")