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
library(patchwork)
library(gridExtra)
library(grid)

# Set seed
set.seed(15)

# Mammal data --------------------------------
# Read in data and sites with no captures
mamm.raw <- read.csv("MammRaw2019.csv")
site.nocaps <- read.csv("AugmentedSites.csv")

# Select columns of interest
mamm.filtered <- mamm.raw %>%
  dplyr::select(Site, Day, Abbrev) %>%
  mutate(Occ = 1)

nocaps.filtered <- site.nocaps %>%
  dplyr::select(Site, Day, Abbrev) %>%
  mutate(Occ = 0)

# Add a species to sites with no caps
mamm.allsite <- rbind(mamm.filtered, nocaps.filtered)
mamm.allsite$Abbrev[is.na(mamm.allsite$Abbrev)] <- "PELE"

# Summary statistics -------------------------
# Number of individuals
inds <- mamm.raw %>%
  dplyr::select(Site, Tag, Abbrev) %>%
  distinct()

nrow(inds)

# Number of species
length(unique(inds$Abbrev))

# Captures by species 
inds %>%
  group_by(Abbrev) %>%
  summarise(byspec = n())

# Prepare data for MSOM input --------------------
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

# Change to presence-absence data
mamm.array[mamm.array > 1] <- 1

# Create array for augmented species
undetected <- array(0, dim = c(J, max(K), 1))

# Put arrays together
mamm.aug <- abind(mamm.array, undetected, along = 3)

# Vegetation data ------------------------
veg <- read.csv("VegRawData.csv")

veg %>%
  dplyr::select(c(Site, Habitat, Canopy:X.DeadVeg)) %>%
  group_by(Site, Habitat) %>%
  summarise_all(mean) %>%
  {. ->> filtered.veg}

# Run PCA
vegpca <- prcomp(filtered.veg[3:15])
#PC1: - grass, + deadveg
#PC2: - forb/bare, + grass
#PC3: - forb, + bare

vegdat <- data.frame(Site = filtered.veg$Site, 
                     Habitat = filtered.veg$Habitat, 
                     PC1 = vegpca$x[,1], PC2 = vegpca$x[,2],
                     PC3 = vegpca$x[,3])

# Plot
pc <- ggplot(data = vegdat, aes(x = scale(PC1), y = scale(PC2), 
                          color = Habitat))+
  geom_point(size = 3)+
  scale_color_viridis_d()+
  labs(x = "PC1", y = "PC2")+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank())+
  annotation_custom(grob = textGrob(label = "Leaf Litter"), 
                    ymin = -2.4, ymax = -2.4, xmin = 1.4, 
                    xmax = 1.4)+
  annotation_custom(grob = textGrob(label = "Grass"), 
                    ymin = -2.4, ymax = -2.4, xmin = -0.8, 
                    xmax = -0.8)+
  annotation_custom(grob = textGrob(label = "Forb/Bare", 
                                    rot = 90),
                    ymin = 1.8, ymax = 1.8, xmin = -1.6, 
                    xmax = -1.6)+
  annotation_custom(grob = textGrob(label = "Grass", rot = 90),
                    ymin = -1.1, ymax = -1.1, xmin = -1.6, 
                    xmax = -1.6)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(pc))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

ggsave(gt, filename = "PC12.jpeg")

# PC1 (separates forest & everything else) explains 82.4% 
forests <- vegdat$PC1
# Standardize
forests <- as.vector(scale(forests))

# PC2 (separates forby & grassy farms/fields) explains 10.5%
farmfield <- vegdat$PC2
# Standardize
farmfield <- as.vector(scale(farmfield))

# Write model script ------------------------------
# Priors
uninf <- "#Add info for species-level priors
            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              w[i] ~ dbern(omega)
              
              a0[i] ~ dnorm(a0.mean, tau.a0)
                             
              a1[i] ~ dnorm(a1.mean, tau.a1)

              b0[i] ~ dnorm(b0.mean, tau.b0)"

weakinf <- "#Add info for species-level priors
            
            inf.mean0 <- -1.3
            inf.mean1 <- -2
            
            inf.var <- 0.5
            
            weights <- c(0.85, 0.15)
            
            lb0[1] <- weights[1]/(1/tau.a0)
            lb0[2] <- weights[2]/inf.var
            lb1[1] <- weights[1]/(1/tau.a1)
            lb1[2] <- weights[2]/inf.var
            
            pooled.var0 <- 1/sum(lb0)
            pooled.mean0 <- sum(lb0*c(a0.mean, inf.mean0))*pooled.var0
            
            pooled.var1 <- 1/sum(lb1)
            pooled.mean1 <- sum(lb1*c(a1.mean, inf.mean1))*pooled.var1
            
            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              w[i] ~ dbern(omega)
              
              a0[i] ~ dnorm(ifelse(i==11, pooled.mean0, a0.mean),
                            ifelse(i==11, (1/pooled.var0), tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==11, pooled.mean1, a1.mean), 
                            ifelse(i==11, (1/pooled.var1), tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

modinf <- "#Add info for species-level priors
            
            inf.mean0 <- -1.3
            inf.mean1 <- -2
            
            inf.var <- 0.5
            
            weights <- c(0.5, 0.5)
            
            lb0[1] <- weights[1]/(1/tau.a0)
            lb0[2] <- weights[2]/inf.var
            lb1[1] <- weights[1]/(1/tau.a1)
            lb1[2] <- weights[2]/inf.var
            
            pooled.var0 <- 1/sum(lb0)
            pooled.mean0 <- sum(lb0*c(a0.mean, inf.mean0))*pooled.var0
            
            pooled.var1 <- 1/sum(lb1)
            pooled.mean1 <- sum(lb1*c(a1.mean, inf.mean1))*pooled.var1
            
            for(i in 1:(spec+aug)){
              #Create priors from hyperpriors
              w[i] ~ dbern(omega)
              
              a0[i] ~ dnorm(ifelse(i==11, pooled.mean0, a0.mean),
                            ifelse(i==11, (1/pooled.var0), tau.a0))
                             
              a1[i] ~ dnorm(ifelse(i==11, pooled.mean1, a1.mean), 
                            ifelse(i==11, (1/pooled.var1), tau.a1))

              b0[i] ~ dnorm(b0.mean, tau.b0)"

# Model text
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
        Z[j,i] ~ dbern(psi[j,i]*w[i])
    
        #Estimate detection of i at point j during sampling period k
        for(k in 1:K[j]){
          logit(p[j,k,i]) <-  b0[i]
          obs[j,k,i] ~ dbern(p[j,k,i]*Z[j,i])
    }
    }
    }
    
    #Estimate total richness by adding observed and unobserved species
    n0<-sum(w[(spec+1):(spec+aug)])
    N<-spec+n0
    
    }
    ")
  writeLines(mod, "realdat.txt") 
}

# Write non-augmented model for testing purposes
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
    
    #Add info for species-level priors
            for(i in 1:spec){
              #Create priors from hyperpriors
              a0[i] ~ dnorm(a0.mean, tau.a0)
                             
              a1[i] ~ dnorm(a1.mean, tau.a1)

              b0[i] ~ dnorm(b0.mean, tau.b0)
    
      #Estimate occupancy of species i at point j
      for (j in 1:J){
        logit(psi[j,i]) <- a0[i] + a1[i]*cov1[j]
        Z[j,i] ~ dbern(psi[j,i])
    
        #Estimate detection of i at point j during sampling period k
        for(k in 1:K[j]){
          logit(p[j,k,i]) <-  b0[i]
          obs[j,k,i] ~ dbern(p[j,k,i]*Z[j,i])
    }
    }
    }
    
    }", file = "noaug.txt")

# Send model to JAGS --------------------------------------------
# Write JAGS function
VivaLaMSOM <- function(J, K, obs, spec = nspec, aug = 1, priors = NULL,
                       cov1 = forests, textdoc = "realdat.txt", 
                       burn = 5000, iter = 15000, thin = 10){
  
  # write the model file
  write.model(priors)
  
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
  parms <- c('N', 'a0.mean', 'a1.mean', 'b0.mean', 'a0', 'b0', 'a1', 
             'Z')
  
  #JAGS command
  model <- jags(model.file = textdoc, data = datalist, 
                n.chains = 3, parameters.to.save = parms, 
                n.burnin = burn, inits = init.values,
                n.iter = iter, n.thin = thin)
  
  return(model)
}

# Test model with no augmented species
# data <- list(J = J, K = K, obs = mamm.array, spec = nspec,
#              cov1 = forests)
# parms <- c('a0.mean', 'a1.mean', 'b0.mean', 'a0', 'b0', 'a1')
# init.values <- function(){
#   inits <- list(Z = apply(mamm.array, c(1,3), max))
# }
# noaug <- jags(model.file = 'noaug.txt', data = data, n.chains = 3,
#               parameters.to.save = parms, n.burnin = 5000,
#               n.iter = 15000, n.thin = 10, inits = init.values)

# Run models and save outputs
# uninf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug, priors = uninf)
# saveRDS(uninf.mod, "real_uninf.rds")
# 
# weakinf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug,
#                           priors = weakinf)
# saveRDS(weakinf.mod, "real_weakinf.rds")
# 
# modinf.mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug, thin = 10,
#                           priors = modinf)
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
  geom_jitter(data = traps, aes(x = X, y = Y, fill = Habitat), 
              size = 3, width = 0.05, height = 0.05, alpha = 0.6,
              pch = 21, color = "black")+
  scale_fill_viridis_d()+
  north(VT, location = "bottomright", scale = 0.2)+
  theme_classic(base_size = 14)+
  theme(axis.title = element_blank(), axis.text = element_blank())

ggsave(file = "sitemap.jpeg", dpi = 600)

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
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), 
        legend.title = element_blank())

ggsave(rich.plot, filename = "SiteRichness_Real.jpeg")

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
  
  full.names <- dplyr::select(mamm.raw, Genus:Abbrev)
  full.names <- rbind(full.names, c("Sylvilagus", "floridanus", "SYFL"))
  
  a1.stat <- a1.long %>%
    group_by(Spec) %>%
    summarise(mean = mean(a1), lo = quantile(a1, 0.025), 
              hi = quantile(a1, 0.975)) %>%
    left_join(y = full.names, by = c("Spec" = "Abbrev")) %>%
    distinct() %>%
    mutate(full.name = paste(substring(Genus,1,1), Species, sep = "."))
  
  # Make interval plot
  plot <- ggplot(data = a1.stat, aes(x = full.name, y = mean))+
    geom_point()+
    geom_errorbar(ymin = a1.stat$lo, ymax = a1.stat$hi, 
                  size = 1, width = 0.2)+
    geom_hline(yintercept = 0, linetype = "dashed", size = 1)+
    scale_y_continuous(limits = c(-12,12))+
    labs(x = "Species", y = "Coefficient")+
    theme_bw(base_size = 14)+
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title = element_blank())
  
  rows <- logical()
  for(i in 1:nrow(a1.stat)){
    rows[i] <- between(0, a1.stat$lo[i], a1.stat$hi[i])
  }
  
  sigs <- rep(NA, nrow(a1.stat))
  sigs[which(rows == FALSE)] <- a1.stat$full.name[which(rows == FALSE)]
  
  outs <- plot+geom_point(aes(x = sigs, y = 6), shape = 8)
  
  return(outs)
}

covslist <- lapply(modlist, get.cov)

all.covplot <- covslist[[1]]/covslist[[2]]/covslist[[3]]+
  plot_annotation(tag_levels = "a")

gcov <- patchworkGrob(all.covplot)
big.covplot <- grid.arrange(gcov, bottom = "Species", left = "Coefficient of Vegetation Cover")

# ggsave(big.covplot, file = "realdatcov.jpeg", width = 6, height = 8,
#        units = 'in')
