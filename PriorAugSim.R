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

# Set seed
set.seed(15)

# Global variables
nspec <- 20
nmiss <- 2 # Species present but not detected during sampling
naug <- 2 # Species never detected; used to set prior on N
nsite <- 30
nsurvey <- 4

Ks <- rep(nsurvey, nsite)

# Priors
uninf <- "for(i in 1:(spec+aug)){
          #Create priors from hyperpriors
            w[i] ~ dbern(omega)
            
            a0[i] ~ dnorm(a0.mean, tau.a0)
            a1[i] ~ dnorm(a1.mean, tau.a1)

            b0[i] ~ dnorm(b0.mean, tau.b0)"


weakinf <- "#Add info for species-level priors
            
            lim <- c(20, 21, 22)

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, -3, 3, 0)
            
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

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
                            
            inf.mean1 <- c(0, -3, 3, 0)
            
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

            inf.mean0 <- c(0, round(logit(sim.occ[21])),
                            round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, -3, 3, 0)
            
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

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, -3, 0)
            
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

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, -3, 0)
            
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

            inf.mean0 <- c(0, -round(logit(sim.occ[21])),
                            -round(logit(sim.occ[22])), 0)
            inf.mean1 <- c(0, 3, -3, 0)
            
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
  resp2cov[21:22] <- c(rnorm(n = 1, mean = -3, sd = 0.25),
                       rnorm(n = 1, mean = 3, sd = 0.25))

  # Covariate values for sites
  cov <- sort(rnorm(n = nsite))
  
  covs.out <- list(resp2cov, cov)
  
  return(covs.out)
}

# Function to simulate occupancy data ---------------------------
occ.func <- function(resp2cov, cov){
  # Get probs from a beta distribution
  sim.occ <- rbeta(n = nspec, shape1 = 2, shape2 = 4)
  # Keep undetected species consistent in all sims
  sim.occ[21:22] <- c(0.3, 0.7)

  #Get site-level psi to account for covariates
  alpha0 <- logit(sim.occ)
  
  logit.psi <- matrix(NA, nrow = nspec+nmiss, ncol = nsite)
  
  for(i in 1:(nspec+nmiss)){
    logit.psi[i,] <- alpha0[i] + resp2cov[i]*cov
  }
  
  psi <- inv.logit(logit.psi)
  
  #create list of abundance vectors
  nlist<-list()
  for(a in 1:(nspec+nmiss)){
    nlist[[a]] <- rbinom(n = nsite, size = 1, prob = psi[a,])
  }
  
  #Turn abundance vectors into abundance matrix
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
  mean.p <- sort(mean.p, decreasing = T)

  # add undetected species
  mean.p[21:22] <- 0

  #Detection intercept and cov responses
  beta0<-logit(mean.p) #put it on logit scale

  #Logit link function
  logit.p <- array(NA, dim = c(nsite, nsurvey, (nspec+nmiss)))
  for(i in 1:(nspec+nmiss)){
    for(j in 1:nsite){
      for(k in 1:nsurvey){
        logit.p[j,,i] <- beta0[i]
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
  
  # Add augmented species
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
                       thin = 10 ){
  
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
  parms <- c('a0', 'a1','Z','N')
  
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
  model <- jags(model.file = textdoc, data = datalist, 
                n.chains = 3, parameters.to.save = parms, 
                inits = init.values, n.burnin = burn, 
                n.iter = iter, n.thin = thin)
  
  return(model)
}

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
    n0<-sum(w[(spec+1):(spec+aug)])
    N<-spec+n0
    
    }
    ")
 writeLines(mod, "aug_model.txt") 
}

# Function to fit models ------------------------------------
fit.mods <- function(cov, obs, sim.occ){
  # Aug model just in case
  # mod.noaug <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
  #                         textdoc = 'noaug.txt')

  mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, spec = nspec, 
                          textdoc = 'aug_model.txt', cov = cov,
                          obs = obs, sim.occ = sim.occ,
                          aug = nmiss+naug, burn = 2500, 
                          iter = 10000, thin = 10)

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
  for(i in 1:iters){
    sim.outs <- comm.sim()
    
    results <- fit.mods(cov = sim.outs[[1]], obs = sim.outs[[3]],
                        sim.occ = sim.outs[[2]])
    
    filename <- paste("./Outputs/out", i, ".rds", sep = "")
    
    saveRDS(results, file = filename)
    
    print(i)
  }
}

# forth.eorlingas(iters = 50)

# Function to load results ---------------
get.outs <- function(param){
  filenames <- list.files("./Outputs")
  
  param.list <- list()
  modtype <- c('mod.uninf', 'inf.weak', 'inf.mod', 'inf.strong',
               'misinf.weak', 'misinf.mod', 'misinf.strong')
  
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
  
  names(param.list) <- modtype
  
  return(param.list)
}

# Compare Ns -------------------------
# Load all n's into R
all.n <- get.outs(param = "N")

# Calculate measure of centrality
ns.center <- function(jag, center = "mean"){
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
            plot.margin = unit(c(0,0,0,0), units = "point"))
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
          plot.margin = unit(c(0,0,0,0), units = "point"))
  
  }
  return(Ns.plot)
}

nouts.mean <- lapply(ns.means, ns.plot, center = "mean")
nouts.mode <- lapply(ns.modes, ns.plot, center = "mode")
nouts.median <- lapply(ns.medians, ns.plot, center = "median")
# I think it's because one spec is rare, another common

# Define plot layout
layout <- "
#AA#
BBEE
CCFF
DDGG
"

histos <- nouts.mode

# Create figure
Ns.base <- histos[[1]]+histos[[2]]+histos[[3]]+histos[[4]]+histos[[5]]+
  histos[[6]]+histos[[7]]+
  plot_layout(design = layout)+
  plot_annotation(tag_levels = "a")+
  plot_layout(guides = "collect")&
  xlim(19.5, 25.5)

gn <- patchworkGrob(Ns.base)
Ns.megaplot <- grid.arrange(gn, bottom = textGrob("Species Richness (N)", 
                                   gp=gpar(fontsize = 14), hjust = 0.8))
  
# ggsave(Ns.megaplot, filename = "ns_median.jpeg", width = 6,
#        height = 5, units = 'in', dpi = 600)

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
  mutate(model = factor(model, 
                        levels = c("mod.uninf", "inf.weak", 
                                   "inf.mod", "inf.strong", 
                                   "misinf.weak", "misinf.mod",
                                   "misinf.strong")))
  
smol.a1.21 <- filter(a1.stat, Species == "Spec21")
smol.a1.22 <- filter(a1.stat, Species == "Spec22")

# Make interval plots
covplot.21 <- ggplot(data = smol.a1.21, aes(x = model, y = mean))+
  geom_point(size = 1.5)+
  geom_errorbar(ymin = smol.a1.21$lo, ymax = smol.a1.21$hi,
                  size = 1, width = 0.2)+
  geom_point(aes(y = -3), color = "red", size = 1.5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0))+
  labs(x = "Model", y = "Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

covplot.22 <- ggplot(data = smol.a1.22, aes(x = model, y = mean))+
  geom_point(size = 1.5)+
  geom_errorbar(ymin = smol.a1.22$lo, ymax = smol.a1.22$hi,
                size = 1, width = 0.2)+
  geom_point(aes(y = 3), color = "red", size = 1.5)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
  scale_y_continuous(limits = c(-10, 10), expand = c(0,0))+
  labs(x = "Model", y = "Coefficient")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())

covs <- (covplot.21/covplot.22)+
  plot_annotation(tag_levels = "a")

# ggsave(covs, filename = "undet_cov.jpeg", dpi = 600, width = 8,
#        height = 6, units = "in")

# Compare site-level richness and covariate ----------------------
# Get true values
sim.res <- list()
for(i in 1:50){
  sim.res[[i]] <- comm.sim()
}

tru <- list()
for(i in 1:length(sim.res)){
  tru[[i]] <- apply(sim.res[[i]][[3]], c(1,3), max)[,-c(23:24)]
}

# Load in estimated vals
zs <- get.outs(param = "Z")

# Get differences between true and est values
get.diff <- function(jag){
  z.avg <- list()
  for(i in 1:7){
    z.avg[[i]] <- list()
    for(j in 1:50){
      z.avg[[i]][[j]] <- apply(jag[[i]][[j]], c(2,3), mean)[,-c(23:24)]
    }
  }

  # Get species richness
  tru.rich <- lapply(tru, rowSums)

  z.rich <- list()
  for(i in 1:length(z.avg)){
      z.rich[[i]] <- lapply(z.avg[[i]], rowSums)
  }

  # subtract true and estimated richness
  z.diff <- list()
  for(i in 1:length(z.rich)){
    z.diff[[i]] <- lapply(z.rich[[i]], function(x) x-tru.rich[[i]])
  }
  names(z.diff) <- names(jag)
  
  return(z.diff)
}

z.diff <- get.diff(jag = zs)

# Get true coefficients
coefs <- list()
for(i in 1: length(sim.res)){
  coefs[[i]] <- sim.res[[i]]$cov
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

modnames <- logical()
for(i in 1:length(z.diff)){
  x <- rep(names(z.diff)[i], 30)
  modnames <- append(modnames, x)
}

# Coerce all to data frames
coef.frame <- as.data.frame(do.call(cbind, coefs))
colnames(coef.frame) <- repnames
coef.frame$Site <- sitenames
coef.frame <- pivot_longer(coef.frame, -Site, names_to = "Rep", 
                           values_to = "Cov")

diff.frame <- lapply(z.diff, function(x) as.data.frame(do.call(cbind, x)))
diff.frame <- do.call(rbind, diff.frame)
diff.frame$Model <- modnames
colnames(diff.frame) <- repnames
diff.frame$Site <- rep(sitenames, 7)
# pivot longer here

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
