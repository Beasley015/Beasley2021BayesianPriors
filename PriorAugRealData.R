###################################################
# Using informed priors with real data            #
# Fall 2019 & Spring 2020                         #
###################################################

# Setup ------------------------------------
# Load packages
library(tidyverse)
library(abind)
library(R2jags)
library(rgdal)
library(spData)
library(sf)
library(ggsn)

# Mammal data --------------------------------
# Read in data and sites with no captures
mamm.raw <- read.csv("MammRaw2019.csv")
site.nocaps <- read.csv("AugmentedSites.csv")

# Select columns of interest
mamm.filtered <- mamm.raw %>%
  select(Site, Day, Abbrev) %>%
  mutate(Occ = 1)

nocaps.filtered <- site.nocaps %>%
  select(Site, Day, Abbrev) %>%
  mutate(Occ = 0)

# Add a species to sites with no caps
mamm.allsite <- rbind(mamm.filtered, nocaps.filtered)
mamm.allsite$Abbrev[is.na(mamm.allsite$Abbrev)] <- "PELE"

# Rearrange data frame
mamm.wide <- mamm.allsite %>%
  complete(Site, Abbrev, Day) %>%
  pivot_wider(names_from = "Day", values_from = "Occ", 
              values_fn = list(Occ = "sum")) 

# Fill missing values with 0
mamm.wide[is.na(mamm.wide)] <- 0

# Get species, sites and vector of sampling periods
specs <- unique(mamm.wide$Abbrev)
nspec <- length(specs)

sites <- unique(mamm.wide$Site)
J <- length(sites)

K <- rep(3, J)

# Coerce mammal data into array
mamm.array <- abind(split(mamm.wide, desc(mamm.wide$Abbrev)), 
                    along = 3)

# Clean it up
mamm.array <- mamm.array[,-(1:2),]
row.names(mamm.array) <- sites
mamm.array <- array(as.numeric(mamm.array), dim = dim(mamm.array))

# Create array for augmented species
undetected <- array(0, dim = c(J, max(K), 1))

# Put arrays together
mamm.aug <- abind(mamm.array, undetected, along = 3)

# Change to presence/absence data
mamm.aug[mamm.aug > 1] <- 1

# Vegetation data ------------------------
veg <- read.csv("VegRawData.csv")

veg %>%
  select(c(Site, Habitat, Canopy:X.DeadVeg)) %>%
  group_by(Site, Habitat) %>%
  summarise_all(mean) %>%
  {. ->> filtered.veg}

# Run PCA
vegpca <- prcomp(filtered.veg[3:15])
#PC1: + grass/forb, - deadveg
#PC2: + forb/bare, - grass/deadveg
#PC3: + forb, - bare

vegdat <- data.frame(Site = filtered.veg$Site, 
                     Habitat = filtered.veg$Habitat, 
                     PC1 = vegpca$x[,1], PC2 = vegpca$x[,2],
                     PC3 = vegpca$x[,3])

# Plot
ggplot(data = vegdat, aes(x = PC1, y = PC2, color = Habitat))+
  geom_point(size = 3)+
  scale_color_viridis_d()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# ggsave("PC12.jpeg")

ggplot(data = vegdat, aes(x = PC2, y = PC3, color = Habitat))+
  geom_point(size = 3)+
  scale_color_viridis_d()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# PC1 (separates forest & everything else) explains 47.5% 
forests <- vegdat$PC1
# Standardize
forests <- as.vector(scale(forests))

# PC2 (separates forby & grassy farms/fields) explains 16%
farmfield <- vegdat$PC2
# Standardize
farmfield <- as.vector(scale(farmfield))

# Write base model script ------------------------------
# Priors
weakinf <- "#Add info for species-level priors
            
            inf.mean0 <- -1.4
            inf.mean1 <- -1
            
            weights <- c(0.85, 0.15)
            
            lb1[1] <- weights[1]/(1/tau.a1)
            lb1[2] <- weights[2]/(1/tau.a1)
              
            pooled.var1 <- 1/sum(lb1)
            
            pooled.mean1 <- sum(lb1*c(a1.mean,inf.mean1))
                               *pooled.var1
            
            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              w[i] ~ dbern(omega)
              
              a0[i] ~ dnorm(a0.mean, tau.a0)
                             
              a1[i] ~ dnorm(ifelse(i==11, pooled.mean1, a1.mean), 
                            ifelse(i==11, (1/pooled.var1), tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

modinf <- "#Add info for species-level priors
            
            inf.mean1 <- -1
            
            inf.var1 <- 0.5
            
            weights <- c(0.5, 0.5)
              
            lb1[1] <- weights[1]/(1/tau.a1)
            lb1[2] <- weights[2]/(1/tau.a1)
              
            pooled.var1 <- 1/sum(lb1)
            
            pooled.mean1 <- sum(lb1*c(a1.mean,inf.mean1))
                               *pooled.var1
            
            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              w[i] ~ dbern(omega)
              
              a0[i] ~ dnorm(a0.mean, tau.a0)
                             
              a1[i] ~ dnorm(ifelse(i==11, pooled.mean1, a1.mean), 
                            ifelse(i==11, (1/pooled.var1), tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

# Model text
uninf.model <- function(){
  mod <- "
    model{
      
    # Define hyperprior distributions: intercepts
    omega ~ dunif(0,1)
    
    #Intercepts
    mean.a0 ~ dunif(0,1)
    a0.mean <- log(mean.a0)-log(1-mean.a0)
    tau.a0 ~ dgamma(0.1, 0.1)
    
    mean.a1 ~ dunif(0,1)
    a1.mean <- log(mean.a1)-log(1-mean.a1)
    tau.a1 ~ dgamma(0.1, 0.1)
    
    mean.b0 ~ dunif(0,1)
    b0.mean <- log(mean.b0)-log(1-mean.b0)
    tau.b0 ~ dgamma(0.1, 0.1)
    
    for(i in 1:(spec+aug)){
      w[i] ~ dbern(omega)
      
      a0[i] ~ dnorm(a0.mean, tau.a0)
      a1[i] ~ dnorm(a1.mean, tau.a1)

      b0[i] ~ dnorm(b0.mean, tau.b0)
    
      #Estimate occupancy of species i at point j
      for (j in 1:J) {
        logit(psi[j,i]) <- a0[i] + a1[i]*cov1[j]
        Z[j,i] ~ dbern(psi[j,i]*w[i])
    
        #Estimate detection of i at point j during sampling period k
        for(k in 1:K[j]){
          logit(p[j,k,i]) <-  b0[i]
          obs[j,k,i] ~ dbern(p[j,k,i]*Z[j,i])
          #The addition of Z means that detecting a species depends on its occupancy
    }
    }
    }
    
    #Estimate total richness (N) by adding observed (n) and unobserved (n0) species
    n0<-sum(w[(spec+1):(spec+aug)])
    N<-spec+n0
    
    }
    "
  writeLines(mod, "realdat_uninf.txt") 
}

# Informed model
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
      for (j in 1:J){
        logit(psi[j,i]) <- a0[i] + a1[i]*cov1[j]
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
  writeLines(mod, "realdat_inf.txt") 
}

# Send model to JAGS --------------------------------------------
# Write JAGS function
VivaLaMSOM <- function(J, K, obs, spec = nspec, aug = 1, priors = NULL,
                       cov1 = forests, textdoc, burn = 2000, 
                       iter = 6000, thin = 5){
  
  # write the model file
  if(textdoc == "realdat_uninf.txt"){
    uninf.model()
  } else{
    write.model(priors)
  }
  
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, aug = aug,
                   cov1 = forests)
  
  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values <- function(){
    omega.guess <- runif(1,0,1)
    inits <- list(
      w = c(rep(1,spec), rbinom(n = aug, size=1, prob=omega.guess)),
      Z = maxobs
    )
  }
  
  # Specify parameters
  parms <- c('N', 'a0.mean', 'b0.mean', 'a0', 'b0', 'a1', 'Z')
  
  #JAGS command
  model <- jags(model.file = textdoc, data = datalist, 
                n.chains = 3, parameters.to.save = parms, 
                n.burnin = burn, inits = init.values,
                n.iter = iter, n.thin = thin)
  
  return(model)
}
# Run models and save outputs
# uninf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug,
#                   textdoc = 'realdat_uninf.txt')
# saveRDS(uninf.mod, "real_uninf.rds")
# 
# weakinf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug,
#                           priors = weakinf,textdoc = "realdat_inf.txt")
# saveRDS(weakinf.mod, "real_weakinf.rds")
# 
# modinf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug,
#                           priors = modinf, textdoc = "realdat_inf.txt")
# saveRDS(modinf.mod, "real_modinf.rds")

# Read in models
uninf.mod <- readRDS("real_uninf.rds")
weakinf.mod <- readRDS("real_weakinf.rds")
modinf.mod <- readRDS("real_modinf.rds")

# Put models in list
modlist <- list(uninf.mod, weakinf.mod, modinf.mod)

# Site map --------------------------
traplines <- readOGR(dsn = "UpdatedTracks.kml")

# Convert to sf object
trapsf <- st_as_sf(traplines)

# Convert to midpoints
st_line_midpoints <- function(sf_lines = NULL){
  
  g <- st_geometry(sf_lines)
  
  g_mids <- lapply(g, function(x){
    
    coords <- as.matrix(x)
    
    # this is just a copypaste of View(maptools:::getMidpoints):
    get_mids <- function (coords){
      dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid <- sum(dist)/2
      dist_cum <- c(0, cumsum(dist))
      end_index <- which(dist_cum > dist_mid)[1]
      start_index <- end_index - 1
      start <- coords[start_index, ]
      end <- coords[end_index, ]
      dist_remaining <- dist_mid - dist_cum[start_index]
      mid <- start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    
    mids <- st_point(get_mids(coords))
  })
}
trap.list <- st_line_midpoints(trapsf)

# Build the data frame
trap.df <- as.data.frame(do.call(rbind, trap.list))
names(trap.df) <- c("X", "Y")
trap.df$Name <- trapsf$Name

# Add habitat types
trap.df %>%
  merge(veg[,c("Site", "Habitat")], by.x = "Name", by.y = "Site") %>%
  distinct() %>%
  {. ->> traps}

# Plot it
VT <- us_states[which(us_states$NAME=="Vermont"),]

ggplot(data = VT)+
  geom_sf()+
  geom_jitter(data = traps, aes(x = X, y = Y, color = Habitat), 
              size = 3, width = 0.05, height = 0.05)+
  scale_color_viridis_d()+
  north(VT, location = "bottomright", scale = 0.2)+
  theme_classic(base_size = 16)+
  theme(axis.title = element_blank(), axis.text = element_blank())

ggsave(file = "sitemap.jpeg")

# Get Ns -----------------------
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
    geom_col(width = 1, color = 'lightgray')+
    geom_vline(xintercept = Ns.median, linetype = 'dashed', size = 2)+
    labs(x = "Estimated Species", y = "Frequency")+
    scale_y_continuous(expand = c(0,0))+
    theme_classic(base_size = 18)+
    theme(axis.text.y = element_blank(), 
          axis.title.y = element_blank(),
          legend.key.height = unit(40, units = 'pt'))
  
  out.list <- list(plot = Ns.plot, mode = Ns.mode, mean = Ns.mean,
                   median = Ns.median)
  
  return(out.list)
}

lapply(modlist, get.ns)

# Site-level richness -------------------
# Pull Zs from each item in list  
Zs <- lapply(modlist, function(x) x$BUGSoutput$sims.list$Z)

# Get avg. occurrence matrices  
Zs.mean <- lapply(Zs, apply, c(2,3), mean)

# Get site-level richness  
site.rich <- lapply(Zs.mean, rowSums)

# Convert list to data.frame
rich.frame <- as.data.frame(do.call(cbind, site.rich))
colnames(rich.frame) <- c('uninf', 'inf.weak', 'inf.mod')

rich.frame$Obs <- rowSums(apply(mamm.aug, c(1,3), max))
rich.frame$Cov <- forests

rich.long <- rich.frame %>%
  pivot_longer(uninf:Obs, names_to = 'model', 
               values_to = 'Richness') %>%
  mutate(mod.type = case_when(startsWith(model,"un")~"Uninformed",
                              endsWith(model,"weak")~"Weakly Informed",
                              endsWith(model,"mod")~"Moderately Informed",
                              model == "Obs"~"Observed"))

order <- c("Observed", "True", "Uninformed", "Weakly Informed",
           "Moderately Informed")
rich.long$mod.type <- factor(rich.long$mod.type, levels = order)

rich.plot <- ggplot(data = rich.long, aes(x = Cov, y = Richness,
                                          color = mod.type))+
  geom_point()+
  geom_smooth(aes(fill = mod.type), method = 'lm', alpha = 0.2)+
  labs(x = "PC1")+
  expand_limits(y = 0)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(), 
        legend.title = element_blank())

# Run individual regressions for fun
mod.names <- unique(rich.long$mod.type)

for(i in seq_along(mod.names)){
  x <- subset(rich.long, rich.long$mod.type == mod.names[i])
  print(summary(lm(data = x, Richness~Cov)))
}

# Covariate responses -------------------
get.cov <- function(jag){
  # Extract covariate estimates from jags object
  a1s <- jag$BUGSoutput$sims.list$a1
  
  a1s <- as.data.frame(a1s)
  
  colnames(a1s) <- c(specs, "SYFL")
  
  # Pivot data frame for plotting
  a1.long <- a1s %>%
    pivot_longer(cols = everything(), names_to = "Spec", 
                 values_to = "a1")
  
  a1.stat <- a1.long %>%
    group_by(Spec) %>%
    summarise(mean = mean(a1), lo = quantile(a1, 0.025), 
              hi = quantile(a1, 0.975))
  
  # Make interval plot
  plot <- ggplot(data = a1.stat, aes(x = Spec, y = mean))+
    geom_point(size = 1.5)+
    geom_errorbar(ymin = a1.stat$lo, ymax = a1.stat$hi, 
                  size = 1, width = 0.2)+
    geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
    scale_y_continuous(limits = c(-5,5))+
    labs(x = "Species", y = "Coefficient")+
    theme_bw(base_size = 14)+
    theme(panel.grid = element_blank())
  
  return(plot)
}

lapply(modlist, get.cov)
