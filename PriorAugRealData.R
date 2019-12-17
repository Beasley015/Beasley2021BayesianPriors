###################################################
# Using informed priors with real data            #
# Fall 2019 & Spring 2020                         #
###################################################

# Setup ------------------------------------
# Load packages
library(tidyverse)

# Load data --------------------------------
# Vegetation data
veg <- read.csv("VegRawData.csv")

veg %>%
  select(c(Site, Habitat, Weins10:Weins60)) %>%
  group_by(Site, Habitat) %>%
  summarise_all(mean) %>%
  {. ->> filtered.veg}

vegpca <- prcomp(filtered.veg[,3:8])

vegdat <- data.frame(Site = filtered.veg$Site, Habitat = filtered.veg$Habitat, 
                     PC1 = vegpca$x[,1], PC2 = vegpca$x[,2])

qplot(data = vegdat, x = PC1, y = PC2, color = Habitat)
