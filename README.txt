README for "Ecologically informed priors improve Bayesian model estimates of species 
richness and occupancy for undetected species"

Beasley 2021

Manuscript currently submitted for review in Ecological Applications

Data: Contains raw mammalian capture data collected throughout Vermont, USA, in summer
2019. Also contains vegetation data collected during the same period, a KML file of 
trapline locations. A .csv file called "AugmentedSites" contains information on sites
where no mammals were captured during the sampling period, which is useful for setting
up data for multi-species occupancy models.

ModelOutputs: Contains posterior distributions (i.e. results) of Bayesian Multi-Species
Occupancy Models described in the manuscript and in Appendix S1.

More detailed metadata is included in each folder.

"PriorAugSim.R" contains the script for generating a simulated metacommunity and testing
effects of prior aggregation on model results.

"PriorAugRealData.R" contains code for applying prior aggregation to empirical data of
small mammal communities in Vermont, USA, collected in summer 2019.
