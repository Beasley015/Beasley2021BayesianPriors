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
mamm.array <- abind(split(mamm.wide, mamm.wide$Abbrev), along = 3)

# Clean it up
mamm.array <- mamm.array[,-(1:2),]
row.names(mamm.array) <- sites
mamm.array <- array(as.numeric(mamm.array), dim = dim(mamm.array))

# Create array for augmented species
undetected <- array(0, dim = c(J, max(K), 2))

# Put arrays together
mamm.aug <- abind(mamm.array, undetected, along = 3)

# Vegetation data ------------------------
veg <- read.csv("VegRawData.csv")

veg %>%
  select(c(Site, Habitat, Weins10:Weins60)) %>%
  group_by(Site, Habitat) %>%
  summarise_all(mean) %>%
  {. ->> filtered.veg}

# Run PCA
vegpca <- prcomp(filtered.veg[,3:8])

vegdat <- data.frame(Site = filtered.veg$Site, Habitat = filtered.veg$Habitat, 
                     PC1 = vegpca$x[,1], PC2 = vegpca$x[,2])

# Plot
ggplot(data = vegdat, aes(x = PC1, y = PC2, color = Habitat))+
  geom_point(size = 3)+
  scale_color_viridis_d()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

# PC1 (lack of structure) explains most 86% of variability- we'll keep that
vegcov <- vegdat$PC1
# Standardize veg covariate
vegcov <- as.vector(scale(vegcov))

# Write base model script ------------------------------
# Priors

# Model text
write.model <- function(priors=NULL){
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
    for(i in 1:(spec+aug)){
      a0[i] ~ dnorm(a0.mean, tau.a0)

      a1[i] ~ dnorm(a1.mean, tau.a1)

      b0[i] ~ dnorm(b0.mean, tau.b0)
    
      #Estimate occupancy of species i at point j
      for (j in 1:J) {
        logit(psi[j,i]) <- a0[i] + a1[i]*cov[j]
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
    ")
  writeLines(mod, "realdat_model.txt") 
}

# Send model to JAGS ------------------------------------------------
# Write JAGS function
VivaLaMSOM <- function(J, K, obs, spec = nspec, aug = 2, cov=vegcov, priors = NULL, 
                       burn = 1000, iter = 5000, thin = 3){
  
  # write the model file
  textdoc <- write.model()
  
  # Compile data into list
  datalist <- list(J = J, K = K, obs = obs, spec = spec, aug = aug, cov = vegcov)
  
  # Specify parameters
  parms <- c('N', 'omega','a0.mean', 'b0.mean', 'a0', 'b0', 'a1','Z')
  
  # Initial values
  maxobs <- apply(obs, c(1,3), max)
  init.values<-function(){
      omega.guess <- runif(1,0,1)
      mu.psi.guess <- runif(1, 0.25, 1)
      inits <- list(
      a0 = rnorm(n = spec+aug), 
      a1 = rnorm(n = spec+aug),
      b0 = rnorm(n = spec+aug),
      Z = maxobs
    )
  }
  
  #JAGS command
  model <- jags(model.file = 'realdat_model.txt', data = datalist, n.chains = 3,
                parameters.to.save = parms, inits = init.values, n.burnin = burn,
                n.iter = iter, n.thin = thin)
  
  return(model)
}

mod <- VivaLaMSOM(J = J, K = K, obs = mamm.aug, cov = vegcov)

# Shapefile --------------------------
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
