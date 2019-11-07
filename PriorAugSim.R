######################################################
# Simulating data for augmentation & informed priors #
# Fall 2019                                          #
######################################################

# Setup -----------------------------------------------
# Load packages
library(R2jags)
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
resp.neutral <- rnorm(n = nspec + naug, mean = 0, sd = 0.25)
resp.present <- c(rnorm(n = 9, mean = 1, sd = 0.25), 
                 rnorm(n = 4, mean = 0, sd = 0.25),
                 rnorm(n = 3, mean = -1, sd = 0.25))

resp2cov <- matrix(c(resp.neutral, resp.present, resp.present), ncol = 3)
resp.names <- c("neutral", "present", "present")

# Covariate values for sites
no.cov <- rep(0, nsite)
cov.bin <- as.numeric(rbernoulli(n = nsite, 0.75))
cov.cont <- rnorm(n = nsite, mean = 0, sd = 2)

cov.vals <- cbind(no.cov, cov.bin, cov.cont)
cov.names <- c("none", "bin", "cont")

# Function for abundance data -------------------------------------
tru.mats <- function(spec = nspec + naug, site = nsite, cov, resp){
  #Draw lambdas from a logseries distribution
  mean.lambdas <- rlogseries(spec, 0.75)
  
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
  
  return(ns)
}

#Get true abundances
for(i in 1:ncol(cov.vals)){
  for(j in 1:ncol(resp2cov)){
     assign(paste(cov.names[i], resp.names[j], sep = ""),
            tru.mats(cov = cov.vals[,i], resp = resp2cov[,j]))
  }
}

trus <- list(binneutral, binpresent, contneutral, contpresent, noneneutral, nonepresent)

#name objects in trus to keep them straight
sim.names <- c('binneutral', 'binpresent', 'contneutral', 'contpresent', 'noneneutral', 
               'nonepresent')

names(trus) <- sim.names

# Simulate detection process ----------------------------------
# Load model results from Master's work
masters.mod <- readRDS("modelsampledglades.rds")

masters.v <- masters.mod$sims.list$v

spec.p <- plogis(masters.v[,1:8])

spec.det <- colMeans(spec.p)

# Get maximum likelihood estimates for params of beta distribution
bet <- fitdistr(x = spec.det, start = list(shape1 = 1, shape2 = 1), "beta")

# Generate detection probabilities from beta dist with above params
sim.dets <- rbeta(n = nspec, shape1 = bet$estimate[1], shape2 = bet$estimate[2])

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

obs.data <- lapply(trus, trap.hist, det = c(sim.dets, 0))

names(obs.data) <- sim.names

# Function to send that sucker to JAGS ----------------------------
# Model without covariate
cat("
    model{
      # Define hyperprior distributions: intercepts
        
        omega ~ dunif(0,1)

        #Intercepts
        a0.mean ~ dnorm(0,0.001)
        sigma.a0 ~ dunif(0,10)
        tau.a0 <- 1/(sigma.a0*sigma.a0)
    
        b0.mean ~ dnorm(0,0.001)
        sigma.b0 ~ dunif(0,10)
        tau.b0 <- 1/(sigma.b0*sigma.b0)
    
        for(i in 1:nspec+naug){
          #create priors from distributions above
          w[i] ~ dbern(omega)
          #indicates whether or not species is exposed to sampling
          
          a0[i] ~ dnorm(a0.mean, tau.a0)
    
          b0[i] ~ dnorm(b0.mean, tau.b0)

        #Estimate occupancy of species i at point j
        for (j in 1:J) {
          logit(psi[j,i]) <- a0[i]+a1[i]*cov[j]
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
        n0<-sum(w[(n+1):(n+zeroes)])
        N<-n+n0

        #Create a loop to determine point level richness estimates
        for(j in 1:J){
          Nsite[j]<- inprod(Z[j,1:(n+zeroes)],w[1:(n+zeroes)])
        }
    }
    ", file = "nocov.txt")
