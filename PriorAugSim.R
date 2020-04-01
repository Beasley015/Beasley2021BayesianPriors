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
library(AICcmodavg)
library(viridis)

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
resp2cov <- c(rnorm(n = 6, sd = 0.25),
              rnorm(n = 5, mean = 2, sd = 0.25),
              rnorm(n = 6, mean = -2, sd = 0.25))

resp2cov <- sample(resp2cov)

detcovresp <- c(rnorm(n = 6, sd = 0.25),
                rnorm(n = 5, mean = 2, sd = 0.25),
                rnorm(n = 6, mean = -2, sd = 0.25))
detcovresp <- sort(detcovresp, decreasing = T)

# Covariate values for sites
cov <- sort(rnorm(n = nsite))

detcov <- rnorm(n = nsite, mean = 1)
detcov <- matrix(rep(detcov, 4), ncol = nsurvey)

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
psi <- rowMeans(tru.mats()[[2]])

# Simulate detection process ---------------------------------
# Generate mean detection probabilities from beta dist
mean.p <- rbeta(n = nspec+nmiss, shape1 = 3, shape2 = 8)
mean.p <- sort(mean.p, decreasing = T)

#There's a problem with this function
get.obs <- function(mat, specs){
  #Detection intercept and cov responses
  beta0<-logit(mean.p) #put it on logit scale

  #Logit link function
  logit.p <- array(NA, dim = c(nsite, nsurvey, specs))
  for(i in 1:specs){
    for(j in 1:nsite){
      for(k in 1:nsurvey){
        logit.p[j,,i] <- beta0[i] + detcovresp[i]*detcov[j,k]
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

#Remove last two specs: these have lowest observations
#And will represent undetected species
obs.data <- obs.data[,,1:nspec]

# Augment the observed dataset ------------------------------------
ems.array <- array(0, dim = c(nsite, nsurvey, ems))
obs.aug <- abind(obs.data, ems.array, along = 3)

# Add prior information --------------------------------
uninf <- list(rep(0, nspec+ems), rep(0, nspec+ems))
weakinf <- list(c(rep(0, nspec), sim.occ[16:17]*0.1, rep(0, naug)),
                  c(rep(0, nspec), round(resp2cov[16:17])*0.1, rep(0, naug)))
modinf <- list(c(rep(0, nspec), sim.occ[16:17]*0.5, rep(0, naug)),
               c(rep(0, nspec), round(resp2cov[16:17])*0.5, rep(0, naug)))
weakmisinf <- list(c(rep(0, nspec), sim.occ[16:17]*-0.1, rep(0, naug)),
                   c(rep(0, nspec), round(resp2cov[16:17])*-0.1, rep(0, naug)))
modmisinf <- list(c(rep(0, nspec), sim.occ[16:17]*-0.5, rep(0, naug)),
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

    mean.b1 ~ dunif(0,1)
    b1.mean <- log(mean.b1)-log(1-mean.b1)
    tau.b1 ~ dgamma(0.1, 0.1)
    
    for(i in 1:spec){
    #create priors from distributions above
    
    a0[i] ~ dnorm(a0.mean, tau.a0)
    a1[i] ~ dnorm(a1.mean, tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    b1[i] ~ dnorm(b1.mean, tau.b1)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    #Estimate detection of i at point j during sampling period k
    for(k in 1:K[j]){
    logit(p[j,k,i]) <-  b0[i] + b1[i]*detcov[j,k]
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

    mean.b1 ~ dunif(0,1)
    b1.mean <- log(mean.b1)-log(1-mean.b1)
    tau.b1 ~ dgamma(0.1, 0.1)
    
    for(i in 1:(spec+aug)){
    # Create priors from hyperpriors
    w[i] ~ dbern(omega)
    #indicates whether or not species is exposed to sampling
    
    a0[i] ~ dnorm(a0.mean+info1[i], tau.a0)
    a1[i] ~ dnorm(a1.mean+info2[i], tau.a1)
    
    b0[i] ~ dnorm(b0.mean, tau.b0)
    b1[i] ~ dnorm(b1.mean, tau.b1)
    
    #Estimate occupancy of species i at point j
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
    mu.psi[j,i] <- psi[j,i] * w[i]
    Z[j,i] ~ dbern(mu.psi[j,i])
    
    #Estimate detection of i at point j during sampling period k
    for(k in 1:K[j]){
    logit(p[j,k,i]) <-  b0[i] + b1[i]*detcov[j,k]
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
  datalist <- list(J = J, K = K, obs = obs, spec = spec, cov = cov, detcov = detcov)
  if(textdoc == 'aug_model.txt'){
    datalist$aug <- aug
  }
  
  if(is.null(info1) == F){
    datalist$info1 <- info1
    datalist$info2 <- info2
  }

  # Specify parameters
  parms <- c('N', 'a0.mean', 'b0.mean', 'a0', 'b0', 'a1', 'b1','Z')

  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
    omega.guess <- runif(1,0,1)
    mu.psi.guess <- runif(1, 0.25, 1)
    inits <- list(
         a0 = rnorm(n = (spec+aug)), a1 = rnorm(n = (spec+aug)),
         b0 = rnorm(n = (spec+aug)), b1 = rnorm(n = (spec+aug)),
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
#            textdoc = 'aug_model.txt', aug = nmiss+naug, info1 = uninf[[1]],
#            info2 = uninf[[2]], burn = 2500, iter = 10000, thin = 10)
# saveRDS(mod.uninf, file = "mod_uninf.rds")
# 
# mod.inf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#                       textdoc = 'aug_model.txt', aug = nmiss+naug,
#                       info1 = weakinf[[1]], info2 = weakinf[[2]], burn = 5000,
#                       iter = 12000, thin = 5)
# saveRDS(mod.inf.weak, file = "mod_inf_weak.rds")
# 
# mod.inf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#                       textdoc = 'aug_model.txt', aug = nmiss+naug,
#                       info1 = modinf[[1]], info2 = modinf[[2]], burn = 7000,
#                       iter = 12000, thin = 3)
# saveRDS(mod.inf, file = "mod_inf.rds")
# 
# mod.misinf.weak <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov,
#                               spec = nspec, textdoc = 'aug_model.txt',
#                               aug = nmiss+naug, info1 = weakmisinf[[1]],
#                               info2 = weakmisinf[[2]], burn = 5000, iter = 10000,
#                               thin = 5)
# saveRDS(mod.misinf.weak, file = "mod_misinf_weak.rds")
# 
# mod.misinf <- VivaLaMSOM(J = nsite, K = Ks, obs = obs.aug, cov = cov, spec = nspec,
#                          textdoc = 'aug_model.txt', aug = nmiss+naug,
#                          info1 = modmisinf[[1]], info2 = modmisinf[[2]], burn = 2500,
#                          iter = 10000, thin = 10)
# saveRDS(mod.misinf, file = "mod_misinf.rds")

# Load models -------------------------------
mod.noaug <- readRDS(file = "mod_noaug.rds")

mod.uninf <- readRDS("mod_uninf.rds")
mod.inf.weak <- readRDS("mod_inf_weak.rds")
mod.inf <- readRDS("mod_inf.rds")
mod.misinf.weak <- readRDS("mod_misinf_weak.rds")
mod.misinf <- readRDS("mod_misinf.rds")

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

  Ns.plot <- ggplot(data = ns.frame, aes(x = as.integer(as.character(N_Species)), 
                                         y = Freq))+
    geom_col()+
    geom_vline(aes(xintercept = Ns.mean, linetype = "Estimated"), size = 1.5)+
    geom_vline(aes(xintercept = nspec+nmiss, linetype = "True"), size = 1.5)+
    scale_linetype_manual(values = c("Estimated"="dotted", "True"="solid"),
                          name = "", labels = c("Mean Estimate", "True"))+
    labs(x = "Number of Species")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_classic(base_size = 18)+
    theme(axis.text.y = element_blank(), legend.key.height = unit(40, units = 'pt'),
          legend.position = c(0.8, 0.8))

  out.list <- list(plot = Ns.plot, mode = Ns.mode, mean = Ns.mean)
  
  return(out.list)
}

N.outs <- lapply(mod.outputs, get.ns)

#View plots
map(N.outs, 1)
map(N.outs, 3)

# Function to compare mean occupancy & detection probabilities ----------------
# Occupancy function
compare.psi <- function(jag){
  psi <- plogis(jag$BUGSoutput$sims.list$a0)
  mean.psi <- plogis(jag$BUGSoutput$sims.list$a0.mean)

  psimat <- data.frame(Observed.Mean = apply(psi, 2, mean)[1:17], 
                      Observed.Lo = apply(psi, 2, quantile, 0.025)[1:17], 
                      Observed.Hi = apply(psi, 2, quantile, 0.975)[1:17],
                      Tru = sim.occ)

  accur <- psimat$Tru >= psimat$Observed.Lo & psimat$Tru <= psimat$Observed.Hi

  perc.acc <- sum(accur)/(nspec+naug)

  psiplot <- ggplot(data = psimat, aes(x = factor(1:17), y = Tru))+
    geom_point(aes(y = Observed.Mean, color = "Estimated"), size = 2)+
    geom_errorbar(aes(ymin = Observed.Lo, ymax = Observed.Hi,
                      color = "Estimated"), size = 1.25)+
    geom_point(aes(color = "True"), size = 2)+
    geom_hline(aes(yintercept = mean(mean.psi)), alpha = 0.5, linetype = 'dashed')+
    scale_color_manual(values = c("black", "red"))+
    labs(x = "Species", y = "Occupancy Probability")+
    theme_bw(base_size = 18)+
    theme(legend.title = element_blank(), panel.grid = element_blank())
  
  outs <- list(accuracy = perc.acc, plot = psiplot)
  
  return(outs)
}

# Detection function
compare.p <- function(jag){
  p <- plogis(jag$BUGSoutput$sims.list$b0)
  mean.p <- plogis(jag$BUGSoutput$sims.list$b0.mean)
  
  pmat <- data.frame(Observed.Mean = apply(p, 2, mean)[1:17], 
                       Observed.Lo = apply(p, 2, quantile, 0.025)[1:17], 
                       Observed.Hi = apply(p, 2, quantile, 0.975)[1:17],
                       Tru = sim.dets)
  
  accur <- pmat$Tru >= pmat$Observed.Lo &
    pmat$Tru <= pmat$Observed.Hi
  
  perc.acc <- sum(accur)/(nspec+naug)
  
  pplot <- ggplot(data = pmat, aes(x = factor(1:17), y = Tru))+
    geom_point(aes(y = Observed.Mean, color = "Estimated"), size = 2)+
    geom_errorbar(aes(ymin = Observed.Lo, ymax = Observed.Hi,
                      color = "Estimated"), size = 1.25)+
    geom_point(aes(color = "True"), size = 2)+
    geom_hline(aes(yintercept = mean(mean.p)), alpha = 0.5, linetype = 'dashed')+
    scale_color_manual(values = c("black", "red"))+
    labs(x = "Species", y = "Detection Probability")+
    theme_bw(base_size = 18)+
    theme(legend.title = element_blank(), panel.grid = element_blank())
  
  outs <- list(accuracy = perc.acc, plot = pplot)
  
  return(outs)
}

# Get outputs
psi.outs <- lapply(mod.outputs, compare.psi)
p.outs <- lapply(mod.outputs, compare.p)

# View accuracy and plots
map(psi.outs, 1)
map(psi.outs, 2)

map(p.outs, 1)
map(p.outs, 2)

# Function to compare covariate responses ----------------------
compare.cov <- function(jag){
  cov.est <- jag$BUGSoutput$sims.list$a1

  covmat <- data.frame(Observed.Mean = apply(cov.est, 2, mean)[1:17], 
                     Observed.Lo = apply(cov.est, 2, quantile, 0.025)[1:17],
                     Observed.Hi = apply(cov.est, 2, quantile, 0.975)[1:17],
                     Tru = resp2cov)

  covplot <- ggplot(data = covmat, aes(x = factor(1:17), y = Tru))+
    geom_point(aes(y = Observed.Mean, color = "Estimated"), size = 2)+
    geom_errorbar(aes(ymin = Observed.Lo, ymax = Observed.Hi,
                      color = "Estimated"), size = 1.25)+
    geom_point(aes(color = "True"), size = 2)+
    geom_hline(aes(yintercept = 0), alpha = 0.5, linetype = 'dashed')+
    scale_color_manual(values = c("black", "red"))+
    labs(x = "Species", y = "Coefficient")+
    theme_bw(base_size = 18)+
    theme(legend.title = element_blank(), panel.grid = element_blank())

  return(covplot)
}

cov.out <- lapply(mod.outputs, compare.cov)

print(cov.out)

# Function to compare true/estimated/observed species richness ------------
richness.comp <- function(jag){
  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.mean <- apply(Zs, c(2,3), mean)

  site.rich <- rowSums(Zs.mean)
  tru.rich <- colSums(tru)
  obs.rich <- rowSums(apply(obs.data, c(1,3), max))

  rich.all <- data.frame(Rank = rank(tru.rich), Estimated = site.rich, 
                         True = tru.rich,Observed = obs.rich)

  rich.plot <- ggplot(data = rich.all, aes(x = Rank, y = True))+
    geom_point(aes(color = "True"))+
    geom_smooth(aes(color = "True"), method = 'lm')+
    geom_point(aes(y = Estimated, color = "Estimated"))+
    geom_smooth(aes(y = Estimated, color = "Estimated"), method = 'lm')+
    geom_point(aes(y = Observed, color = "Observed"))+
    geom_smooth(aes(y = Observed, color = "Observed"), method = 'lm')+
    labs(x = "Sites (Ranked)")+
    expand_limits(y = 0)+
    scale_color_viridis_d()+
    theme_bw(base_size = 18)+
    theme(panel.grid = element_blank(), legend.title = element_blank())
  
  return(rich.plot)
}

rich.out <- lapply(mod.outputs, richness.comp)

# Function looking at observed~true occupancy -------------------------
error.raster <- function(jag){
  specnames <- as.character(1:(nspec+nmiss))

  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.mean <- data.frame(apply(Zs, c(2,3), mean)[,1:17])
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

# Compare true vs. estimated richness ------------------------
richness.bias <- function(jag){
  Zs <- jag$BUGSoutput$sims.list$Z
  
  Zs.mean <- apply(Zs, c(2,3), mean)
  
  site.rich <- rowSums(Zs.mean)
  tru.rich <- colSums(tru)
  obs.rich <- rowSums(apply(obs.data, c(1,3), max))
  
  rich.all <- data.frame(Estimated = site.rich, True = tru.rich, 
                         True.Sq = tru.rich^2)
  
  #Plot examining relationship
  rich.plot <- ggplot(data = rich.all, aes(x = True, y = Estimated))+
    geom_point()+
    geom_abline(slope = 1, intercept = 0)+
    geom_smooth(color = "black", method = 'lm')+
    theme_bw(base_size = 18)+
    theme(panel.grid = element_blank())
  
  # Examine fits of different curves
  linear <- glm(data = rich.all, Estimated~True)
  expon <- glm(data = rich.all, Estimated~True.Sq)
  quad <- glm(data = rich.all, Estimated~True+True.Sq)
  
  comparison <- aictab(cand.set = list(linear, expon, quad),
                       modnames = c("linear", "exponential", "quadratic"))
  
  outs <- list(Plot = rich.plot, AIC = comparison, Model = summary(linear))
  
  return(outs)
}

bias.out <- lapply(mod.outputs, richness.bias)
map(bias.out, 1)
map(bias.out, 2)
map(bias.out, 3)

# Compare error and detection probability ---------------------------
erdet <- function(jag){
  specnames <- as.character(1:(nspec+nmiss))

  Zs <- jag$BUGSoutput$sims.list$Z

  Zs.mean <- data.frame(apply(Zs, c(2,3), mean)[,1:17])
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
    merged.frame$Detection[i] <- sim.dets[as.numeric(merged.frame$Species[i])]
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

