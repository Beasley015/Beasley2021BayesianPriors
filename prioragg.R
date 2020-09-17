# Logarithmic pooling of priors (Gaussian)
# Fall 2020
# From Leo Bastos & Luiz Max Carvalho (2019) and Malta et al. (2010)

library(tidyverse)

# Define functions
pool_par_gauss <- function(alpha, m, v){
  ## same as pool_par(), but for the Gaussian distribution
  ## alpha = vector of weights
  ## Takes in MEANS and VARIANCES and outputs MEAN and SD (for plotting reasons)
  ws <- alpha/v
  vstar <-  1/sum(ws)
  mstar <- sum(ws*m) * vstar
  c(mstar, sqrt(vstar))
}

stat_gauss <- function(p, alpha = .95){
  # returns the mean and 100*alpha% quantiles
  c( p[1], qnorm( c((1-alpha)/2, (1+alpha)/2), p[1], p[2]))
}

############
K <- 2 # Number of priors to aggregate
mv <- c(0,3) # Vector of means
sv <- c(0.5, 1) # Vector of standard deviations

# Equal weights
alphaEqual <- rep(1/K, K)

ab.Equal.star <- pool_par_gauss(alphaEqual, mv, sv^2)

# Plot
ggplot(mapping = aes(x = c(-2, 5)))+
  stat_function(fun = dnorm, n = 100, 
                args = list(mean = ab.Equal.star[1], sd = ab.Equal.star[2]))+
  stat_function(fun = dnorm, n = 100, args = list(mean = mv[1], sd = sv[1]), 
                linetype = "dashed")+
  stat_function(fun = dnorm, n = 100, args = list(mean = mv[2], sd = sv[2]), 
                linetype = "dashed")+
  labs(y = "Density", x = "X")+
  theme_bw()+
  theme(panel.grid = element_blank())

  

