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

# Set seed
set.seed(16)

# Global variables
nspec <- 20
nmiss <- 2 # Species present but not detected during sampling
naug <- 3 # Hypothetical species not exposed to sampling
ems <- nmiss + naug # Number of all-zero histories to add to observed data
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

#Sanity check: get observed data matrix
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
ems.array <- array(0, dim = c(nsite, nsurvey, ems))
obs.aug <- abind(obs.data, ems.array, along = 3)

# Add prior information --------------------------------
uninf <- "for(i in 1:(spec+aug)){
          #Create priors from hyperpriors
            w[i] ~ dbern(omega)
            #indicates whether or not species is exposed to sampling

            a0[i] ~ dnorm(a0.mean, tau.a0)
            a1[i] ~ dnorm(a1.mean, tau.a1)

            b0[i] ~ dnorm(b0.mean, tau.b0)"


weakinf <- "#Add info for species-level priors
            lim <- c(20, 21, 22)    
            
            #Intercept information
            a0.lo1 <- a0.mean-(5/sqrt(tau.a0)) 
            a0.lo2 <- a0.mean-(2/sqrt(tau.a0))
            a0.lo3 <- a0.mean-(5/sqrt(tau.a0))
              
            a0.lo <- c(a0.lo3, a0.lo1, a0.lo2, a0.lo3)
            
            a0.hi1 <- a0.mean+(1/sqrt(tau.a0))
            a0.hi2 <- a0.mean+(2/sqrt(tau.a0))
            a0.hi3 <- a0.mean+(5/sqrt(tau.a0))
            
            a0.hi <- c(a0.hi3, a0.hi1, a0.hi2, a0.hi3)
            
            #Covariate information
            a1.lo1 <- a1.mean-(2/sqrt(tau.a1))
            a1.lo2 <- a1.mean-(1/sqrt(tau.a1))
            a1.lo3 <- a1.mean-(5/sqrt(tau.a1))
            
            a1.lo <- c(a1.lo3, a1.lo1, a1.lo2, a1.lo3)
            
            a1.hi1 <- a1.mean+(2/sqrt(tau.a1))
            a1.hi2 <- a1.mean+(5/sqrt(tau.a1))
            a1.hi3 <- a1.mean+(5/sqrt(tau.a1))
            
            a1.hi <- c(a1.hi3, a1.hi1, a1.hi2, a1.hi3)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              
              w[i] ~ dbern(ifelse(i == 21 || i == 22, 0.75, omega))
              #indicates whether or not species is exposed to sampling"

# Note: strongly informed priors below need inits close to the true value to work
modinf <- "#Add info for species-level priors
            lim <- c(20, 21, 22)    
            
            #Intercept information
            a0.lo1 <- a0.mean-(5/sqrt(tau.a0)) 
            a0.lo2 <- a0.mean-(1/sqrt(tau.a0))
            a0.lo3 <- a0.mean-(5/sqrt(tau.a0))
              
            a0.lo <- c(a0.lo3, a0.lo1, a0.lo2, a0.lo3)
            
            a0.hi1 <- a0.mean+0.001
            a0.hi2 <- a0.mean+(2/sqrt(tau.a0))
            a0.hi3 <- a0.mean+(5/sqrt(tau.a0))
            
            a0.hi <- c(a0.hi3, a0.hi1, a0.hi2, a0.hi3)
            
            #Covariate information
            a1.lo1 <- a1.mean-(2/sqrt(tau.a1))
            a1.lo2 <- a1.mean-0.001
            a1.lo3 <- a1.mean-(5/sqrt(tau.a1))
            
            a1.lo <- c(a1.lo3, a1.lo1, a1.lo2, a1.lo3)
            
            a1.hi1 <- a1.mean+(2/sqrt(tau.a1))
            a1.hi2 <- a1.mean+(5/sqrt(tau.a1))
            a1.hi3 <- a1.mean+(5/sqrt(tau.a1))
            
            a1.hi <- c(a1.hi3, a1.hi1, a1.hi2, a1.hi3)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              
              w[i] ~ dbern(ifelse(i == 21 || i == 22, 0.95, omega))
              #indicates whether or not species is exposed to sampling"
  
# Inits also need adjusting here
weakmisinf <- "#Add info for species-level priors
            lim <- c(20, 21, 22)    
            
            #Intercept information
            a0.lo1 <- a0.mean-0.001 
            a0.lo2 <- a0.mean-(2/sqrt(tau.a0))
            a0.lo3 <- a0.mean-(5/sqrt(tau.a0))
              
            a0.lo <- c(a0.lo3, a0.lo1, a0.lo2, a0.lo3)
            
            a0.hi1 <- a0.mean+(5/sqrt(tau.a0))
            a0.hi2 <- a0.mean+(1/sqrt(tau.a0))
            a0.hi3 <- a0.mean+(5/sqrt(tau.a0))
            
            a0.hi <- c(a0.hi3, a0.hi1, a0.hi2, a0.hi3)
            
            #Covariate information
            a1.lo1 <- a1.mean-(5/sqrt(tau.a1))
            a1.lo2 <- a1.mean-(5/sqrt(tau.a1))
            a1.lo3 <- a1.mean-(5/sqrt(tau.a1))
            
            a1.lo <- c(a1.lo3, a1.lo1, a1.lo2, a1.lo3)
            
            a1.hi1 <- a1.mean+(1/sqrt(tau.a1))
            a1.hi2 <- a1.mean+0.001
            a1.hi3 <- a1.mean+(5/sqrt(tau.a1))
            
            a1.hi <- c(a1.hi3, a1.hi1, a1.hi2, a1.hi3)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              
              w[i] ~ dbern(ifelse(i == 21 || i == 22, 0.25, omega))
              #indicates whether or not species is exposed to sampling"
 
#Inits adjusted here too 
modmisinf <- "#Add info for species-level priors
            lim <- c(20, 21, 22)    
            
            #Intercept information
            a0.lo1 <- a0.mean+(1/sqrt(tau.a0))
            a0.lo2 <- a0.mean-(1/sqrt(tau.a0))
            a0.lo3 <- a0.mean-(5/sqrt(tau.a0))
              
            a0.lo <- c(a0.lo3, a0.lo1, a0.lo2, a0.lo3)
            
            a0.hi1 <- a0.mean+(5/sqrt(tau.a0))
            a0.hi2 <- a0.mean+(1/sqrt(tau.a0))
            a0.hi3 <- a0.mean+(5/sqrt(tau.a0))
            
            a0.hi <- c(a0.hi3, a0.hi1, a0.hi2, a0.hi3)
            
            #Covariate information
            a1.lo1 <- a1.mean-(5/sqrt(tau.a1))
            a1.lo2 <- a1.mean-(5/sqrt(tau.a1))
            a1.lo3 <- a1.mean-(5/sqrt(tau.a1))
            
            a1.lo <- c(a1.lo3, a1.lo1, a1.lo2, a1.lo3)
            
            a1.hi1 <- a1.mean+0.001
            a1.hi2 <- a1.mean-(1/sqrt(tau.a1))
            a1.hi3 <- a1.mean+(5/sqrt(tau.a1))
            
            a1.hi <- c(a1.hi3, a1.hi1, a1.hi2, a1.hi3)

            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              g[i] ~ dinterval(i, lim)
              
              w[i] ~ dbern(ifelse(i == 21 || i == 22, 0.25, omega))
              #indicates whether or not species is exposed to sampling"

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

      a0[i] ~ dnorm(a0.mean, tau.a0)T(a0.lo[g[i]+1], a0.hi[g[i]+1])

      a1[i] ~ dnorm(a1.mean, tau.a1)T(a1.lo[g[i]+1], a1.hi[g[i]+1])

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
    ")
 writeLines(mod, "aug_model.txt") 
}

# Write function for sending model to gibbs sampler --------------------------------
VivaLaMSOM <- function(J, K, obs, spec, aug = 0, cov, textdoc, priors = uninf, 
                       burn = 2500, iter = 8000, thin = 10){
  # Write model for augmented datasets
  if(textdoc == 'aug_model.txt')
    write.model(priors = priors)
  
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, cov = cov)
  if(textdoc == 'aug_model.txt'){
    datalist$aug <- aug
  }

  # Specify parameters
  parms <- c('N', 'omega','a0.mean', 'b0.mean', 'a0', 'b0', 'a1','Z')

  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
    omega.guess <- runif(1,0,1)
    mu.psi.guess <- runif(1, 0.25, 1)
    inits <- list(
         a0 = c(-logit(sim.occ)+rnorm(nspec+nmiss,0,0.1),rnorm(naug)), 
         a1 = c(-resp2cov+rnorm(nspec+nmiss, 0, 0.1), rnorm(naug)),
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
# 
# mod.uninf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#            textdoc = 'aug_model.txt', aug = nmiss+naug, burn = 2500, iter = 10000,
#            thin = 10)
# saveRDS(mod.uninf, file = "mod_uninf.rds")
# 
# mod.inf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                            spec = nspec,textdoc = 'aug_model.txt',
#                            aug = nmiss+naug, priors = weakinf, burn = 5000,
#                            iter = 12000, thin = 5)
# saveRDS(mod.inf.weak, file = "mod_inf_weak.rds")
# 
# mod.inf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#                       textdoc = 'aug_model.txt', aug = nmiss+naug, priors = modinf,
#                       burn = 7000, iter = 12000, thin = 3)
# saveRDS(mod.inf, file = "mod_inf.rds")
# 
# mod.misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                               spec = nspec, textdoc = 'aug_model.txt',
#                               aug = nmiss+naug, priors = weakmisinf,
#                               burn = 5000, iter = 10000, thin = 5)
# saveRDS(mod.misinf.weak, file = "mod_misinf_weak.rds")
# 
# mod.misinf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#                          textdoc = 'aug_model.txt', aug = nmiss+naug,
#                          priors = modmisinf, burn = 2500, iter = 10000, thin = 10)
# saveRDS(mod.misinf, file = "mod_misinf.rds")

# Load models -------------------------------
mod.noaug <- readRDS(file = "mod_noaug.rds")

# Models with uninformed covariates
mod.uninf <- readRDS("mod_uninf.rds")
mod.inf.weak <- readRDS("wionly/mod_inf_weak.rds")
mod.inf <- readRDS("wionly/mod_inf.rds")
mod.misinf.weak <- readRDS("wionly/mod_misinf_weak.rds")
mod.misinf <- readRDS("wionly/mod_misinf.rds")

# Models with informed covariates
mod.uninf <- readRDS("mod_uninf.rds")
mod.inf.weak <- readRDS("mod_inf_weak.rds")
mod.inf <- readRDS("mod_inf.rds")
mod.misinf.weak <- readRDS("mod_misinf_weak.rds")
mod.misinf <- readRDS("mod_misinf.rds")

# Put models in a list
mod.outputs <- list(mod.uninf, mod.inf.weak, mod.inf, mod.misinf.weak, mod.misinf)

# Function to compare estimates of N -----------------------------------
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

  Ns.plot <- ggplot(data = ns.frame, aes(x = as.integer(as.character(N_Species)), 
                                         y = Freq))+
    geom_col()+
    geom_vline(aes(xintercept = Ns.median, linetype = "Estimated"), size = 1.5)+
    geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"), size = 1.5)+
    scale_linetype_manual(values = c("Estimated"="dotted", "True"="solid"),
                          name = "", labels = c("Median Estimate", "True"))+
    labs(x = "Estimated Species", y = "Frequency")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_classic(base_size = 18)+
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          legend.key.height = unit(40, units = 'pt'))

  out.list <- list(plot = Ns.plot, mode = Ns.mode, mean = Ns.mean, median = Ns.median)
  
  return(out.list)
}

N.outs <- lapply(mod.outputs, get.ns)

#View plots
map(N.outs, 1)
map(N.outs, 2)
map(N.outs, 3)
map(N.outs, 4)

# Put histograms in single figure
layout <- 
  "#AA#
   BBCC
   DDEE"

histos <- map(N.outs, 1) 

Ns.megaplot <- histos[[1]]+ ggtitle('A)') +
  histos[[2]] + ggtitle('B)') +
  histos[[3]] + ggtitle('C)') +
  histos[[4]] + ggtitle('D)') +
  histos[[5]] + ggtitle('E)') +
  plot_layout(design = layout, guides = 'collect')

# Save figure
# ggsave(Ns.megaplot, file = "Nswionly.jpeg", height = 7,width = 7, units = "in")


# Get omegas ------------------------------------------
get.omega <- function(jag){
  os <- jag$BUGSoutput$sims.list$omega
  
  omega.plot <- ggplot()+
    geom_density(aes(x = os), fill = 'lightgray')+
    labs(x = "Regional Occurrence Probability", y = "Density")+
    scale_x_continuous(limits = c(0,1))+
    theme_bw(base_size = 18)+
    theme(panel.grid = element_blank())
  
  return(omega.plot)
}

omegas.out <- lapply(mod.outputs, get.omega)

Os.megaplot <- omegas.out[[1]]+ ggtitle('A)') +
  omegas.out[[2]] + ggtitle('B)') +
  omegas.out[[3]] + ggtitle('C)') +
  omegas.out[[4]] + ggtitle('D)') +
  omegas.out[[5]] + ggtitle('E)') +
  plot_layout(design = layout)

# Compare site-level richness ---------------------------
# Create list with unaugmented data
biglist <- list(mod.noaug, mod.uninf, mod.inf.weak, mod.inf, mod.misinf.weak, 
                mod.misinf)

# Pull Zs from each item in list  
Zs <- lapply(biglist, function(x) x$BUGSoutput$sims.list$Z)

# Get avg. occurrence matrices  
Zs.mean <- lapply(Zs, apply, c(2,3), mean)

# Get site-level richness  
site.rich <- lapply(Zs.mean, rowSums)

# Convert list to data.frame
rich.frame <- as.data.frame(do.call(cbind, site.rich))
colnames(rich.frame) <- c('noaug', 'uninf', 'inf.weak', 'inf', 'misinf.weak',
                          'misinf')

rich.frame$True <- colSums(tru)
rich.frame$Obs <- rowSums(apply(obs.data, c(1,3), max))
rich.frame$Rank <- rank(rich.frame$True, ties.method = 'first')

rich.long <- rich.frame %>%
  pivot_longer(noaug:Obs, names_to = 'model', values_to = 'Richness')
  
rich.plot <- ggplot(data = rich.long, aes(x = Rank, y = Richness, color = model))+
  geom_point()+
  geom_smooth(aes(fill = model), method = 'lm', alpha = 0.2)+
  labs(x = "Sites (Ranked)")+
  expand_limits(y = 0)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank(), legend.title = element_blank())

ggsave(rich.plot, filename = 'richplot.jpeg')

# Compare typical bias of each prior method ---------------------------
# Get series of site-level estimates from Zs
Zs.replicates <- lapply(Zs, apply, c(1,2), sum)

# Subtract each replicate vector from the true value
rich.diffs <- lapply(Zs.replicates, apply, 1, function(x) x-colSums(tru))

# Convert each matrix in list into a giant vector
rich.vecs <- lapply(rich.diffs, as.vector)

# Get avg deviance from true value
mean.diff <- unlist(lapply(rich.diffs, mean))

# Create histograms
make.rich.hists <- function(x){
  plot <- ggplot()+
    geom_bar(aes(x = x), fill = 'darkgray')+
    geom_vline(xintercept = 0)+
    labs(x = "Richness")+
    theme_bw(base_size = 18)+
    theme(axis.title.y = element_blank(), panel.grid = element_blank())
  
  return(plot)
}

rich.hists <- lapply(rich.vecs, make.rich.hists)

layout2 <- "
AABB
CCDD
EEFF"

allthehists <- rich.hists[[1]]+ ggtitle("No Augmentation")+
  rich.hists[[2]]+ ggtitle("Uninformed")+
  rich.hists[[3]]+ ggtitle("Weakly Informed")+
  rich.hists[[4]]+ ggtitle("Informed")+
  rich.hists[[5]]+ ggtitle("Weakly Misinformed")+
  rich.hists[[6]]+ ggtitle("Misinformed")+
  plot_layout(design = layout2)

# ggsave(allthehists, file = "allthehists.jpeg", height = 8, width = 8, units = 'in')

# Function to compare covariate responses ----------------------
a1s <- lapply(biglist, function(x) x$BUGSoutput$sims.list$a1)

a1.frame <- lapply(a1s, as.data.frame)

specnames <- paste("Spec", c(1:25), sep = "")

a1.frame <- lapply(a1.frame, function(x) setNames(x,specnames[1:ncol(x)]))

longggggg.frame <- function(x){
  y <- x %>%
    pivot_longer(cols = everything(), names_to = "Spec", values_to = "a1")
  
  return(y)
}

a1.long <- lapply(a1.frame, longggggg.frame)

trim <- function(x){
  a1.trim <- x %>%
    group_by(Spec) %>%
    filter(between(a1, quantile(a1, 0.025), quantile(a1, 0.975)))

  return(a1.trim)
}

a1.trim <- lapply(a1.long, trim)

make.violins <- function(dat){
  ggplot(data = dat, aes(x = Spec, y = a1))+
    geom_violin(fill = 'lightgray', draw_quantiles = c(0.025, 0.975))+
    #geom_point(aes(y = resp2cov))+
    geom_hline(yintercept = 0, linetype = "dashed", size = 1.5)+
    theme_bw(base_size = 18)+
    theme(axis.title.y = element_blank(), panel.grid = element_blank())
}

lapply(a1.trim, make.violins)

# Function looking at observed~true occupancy -------------------------
error.raster <- function(jag){
  specnames <- as.character(1:(nspec+nmiss))

  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.mean <- data.frame(apply(Zs, c(2,3), mean)[,1:22])
  colnames(Zs.mean) <- specnames
  Zs.mean$Site <- 1:nrow(Zs.mean)

  tru.frame <-as.data.frame(t(tru))
  colnames(tru.frame) <- specnames
  tru.frame$Site <- 1:nrow(tru.frame)

  tru.frame %>%
    gather('1':'17', key = "Species", value = "Occ") %>%
    {. ->> tru.frame}

  Zs.mean %>%
    gather('1':'17', key = "Species", value = "Occ") %>%
    left_join(tru.frame, by = c("Site", "Species")) %>%
    mutate(Error = Occ.y - Occ.x) %>%
    {. ->> merged.frame}

  merged.frame$Species = factor(merged.frame$Species, levels = specnames)

  rast <- ggplot(data = merged.frame, aes(x = Species, y = Site, fill = Error))+
    geom_tile()+
    scale_fill_distiller(type = "div", 
                         limit = max(abs(merged.frame$Error)) * c(-1, 1))+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    theme_bw(base_size = 18)
  
  return(rast)
}

error.out <- lapply(mod.outputs, error.raster)

print(error.out)

# Compare error and detection probability ---------------------------
erdet <- function(jag){
  specnames <- as.character(1:(nspec+nmiss))

  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.mean <- data.frame(apply(Zs, c(2,3), mean)[,1:22])
  colnames(Zs.mean) <- specnames
  Zs.mean$Site <- 1:nrow(Zs.mean)

  tru.frame <-as.data.frame(t(tru))
  colnames(tru.frame) <- specnames
  tru.frame$Site <- 1:nrow(tru.frame)

  tru.frame %>%
    gather('1':'17', key = "Species", value = "Occ") %>%
    {. ->> tru.frame}

  Zs.mean %>%
    gather('1':'17', key = "Species", value = "Occ") %>%
    left_join(tru.frame, by = c("Site", "Species")) %>%
    mutate(Error = Occ.y - Occ.x) %>%
    {. ->> merged.frame}

  for(i in 1:nrow(merged.frame)){
    merged.frame$Detection[i] <- mean.p[as.numeric(merged.frame$Species[i])]
  }

  hexplot <- ggplot(data = merged.frame, aes(x = Detection, y = Error))+
    stat_binhex()+ 
    geom_hline(yintercept = 0)+
    scale_fill_viridis(name = "Count")+
    labs(x = "Detection Probability")+
    theme_bw(base_size = 20)+
    theme(panel.grid = element_blank())
  
  return(hexplot)
}

erdet.out <- lapply(mod.outputs, erdet)

print(erdet.out)
