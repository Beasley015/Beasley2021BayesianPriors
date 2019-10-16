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

