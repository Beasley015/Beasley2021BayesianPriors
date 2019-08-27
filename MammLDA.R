###########################################################
#  Classification of Small Mammal Communities in VT       #
#  MSAM + LDA                                             #
#  Summer/Fall 2019                                       #
###########################################################

# Load packages and set wd -------------------------------------------

library(R2OpenBUGS)
library(tidyverse)
library(reshape2)

# Load and clean mammal data -----------------------------------------
mamm <- read.csv("mammraw2019.csv", stringsAsFactors = F)

# Trim dataset to include species/sites/surveys
mamm.smol <- select(mamm, Site, Day, Abbrev)

# Add sites with 0 captures
no.caps <- data.frame(Site = c("Borderview2", "Buck1", "Buck3", "Buck4", "Butternut1",
                               "MBR2", "MBR3", "MBR4", "RiverBerry2"),
                      Day = rep(1, 9),
                      Abbrev = rep("PELE", 9))

mamm.smol <- rbind(mamm.smol, no.caps)

# Get abundances per site/species/day
mamm.smol %>%
  group_by(Site, Day, Abbrev) %>%
  summarise(Abund = length(Site)) %>%
  {. ->> mamm.clean}

# Convert above to 3D array to send to model
mamm.input <- acast(mamm.clean, formula = Site~Day~Abbrev)
mamm.input[is.na(mamm.input)] <- 0

# Create abundance matrix from dataframe
mamm.clean %>%
  ungroup() %>%
  group_by(Abbrev, Site) %>%
  summarise(MaxAbund = max(Abund)) %>%
  spread(key = Site, value = MaxAbund) %>%
  {. ->> abund.mat}

abund.mat[is.na(abund.mat)] <- 0
abund.mat <- abund.mat[,-1]

# Get number of species/site/samples ------------------------------- 
# Get number of unique species and sites
unique.specs <- as.numeric(length(unique(mamm.clean$Abbrev)))
unique.sites <- as.numeric(length(unique(mamm.clean$Site)))

# Create vector of number of days per site
# There were 3 sampling periods per site so this is an easy one
K <- rep(3, unique.sites)

# Run the MSAM -----------------------------------
cat("
    model{
    
      # Define prior distributions

      # Intercepts
      a0.mean ~ dnorm(0,0.001)
      sigma.a0 ~ dunif(0,10)
      tau.a0 <- 1/(sigma.a0*sigma.a0)
    
      b0.mean ~ dnorm(0,0.001)
      sigma.b0 ~ dunif(0,10)
      tau.b0 <- 1/(sigma.b0*sigma.b0)
      
      # Covariates would go here

      for(i in 1:N){
        # Create priors from distributions above
        a0[i] ~ dnorm(a0.mean, tau.a0)
            
        b0[i] ~ dnorm(b0.mean, tau.b0)

        # Estimate abundance of species i at site j
        for(j in 1:J){
          log(lambda[j,i]) <- a0[i]
          Z[j,i] ~ dpois(lambda[j,i])

          # Estimate detection of spec i at site j during sample k
          for(k in 1:K[j]){
            p[j,k,i] <- b0[i]
            logit.p[j,k,i] <- 1 / (1 + exp(-p[j,k,i]))
            obs[j,k,i] ~ dbin(logit.p[j,k,i], Z[j,i])
        }
      }
    }
  }
    ", file = "MammCommModel.txt")

# Compile data into list
datalist <- list(N = unique.specs, J = unique.sites, K = K, obs = mamm.input)

# Specify parameters to return
params <- list('Z', 'lambda','a0','b0', 'a0.mean', 'b0.mean')

# Specify initial values
init.values <- function(){
  list(a0 = rnorm(n = unique.specs), b0 = rnorm(n = unique.specs),
       Z = t(as.matrix(abund.mat)))
}

# Send model to Gibbs sampler
# mod <- bugs(data = datalist, inits = init.values, parameters.to.save = params,
#             n.iter = 5000, model.file = "MammCommModel.txt", n.chains = 3,
#             n.burnin = 2000, debug = T)
# 
# saveRDS(mod, file = "MSAMoutput2019.rds")

mod <- readRDS(file = "MSAMoutput2019.rds")

# Summary -------------------------------------------------
Zs <- mod$sims.list$Z
Zmean <- as.data.frame(apply(Zs, c(2,3), mean))

colnames(Zmean) <- sort(unique(mamm.clean$Abbrev))

