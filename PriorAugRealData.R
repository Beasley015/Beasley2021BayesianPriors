###################################################
# Using informed priors with real data            #
# Fall 2019 & Spring 2020                         #
###################################################

# Setup ------------------------------------
# Load packages
library(tidyverse)
library(rgdal)
library(spData)
library(sf)
library(ggsn)

# Mammal data --------------------------------
mamm.raw <- read.csv("MammRaw2019.csv")
site.nocaps <- read.csv("AugmentedSites.csv")

mamm.filtered <- mamm.raw %>%
  select(Site, Day, Abbrev)

# Vegetation data ------------------------
veg <- read.csv("VegRawData.csv")

veg %>%
  select(c(Site, Habitat, Weins10:Weins60)) %>%
  group_by(Site, Habitat) %>%
  summarise_all(mean) %>%
  {. ->> filtered.veg}

vegpca <- prcomp(filtered.veg[,3:8])

vegdat <- data.frame(Site = filtered.veg$Site, Habitat = filtered.veg$Habitat, 
                     PC1 = vegpca$x[,1], PC2 = vegpca$x[,2])

ggplot(data = vegdat, aes(x = PC1, y = PC2, color = Habitat))+
  geom_point(size = 3)+
  scale_color_viridis_d()+
  theme_bw(base_size = 18)+
  theme(panel.grid = element_blank())

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
