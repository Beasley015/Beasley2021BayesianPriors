#######################################################
# Plot of hypothetical scenarios                      #
# Fig. 1 of proposal                                  #
# Spring 2020                                         #
#######################################################

# Load package
library(ggplot2)

# Write function for logistic curve
curve <- function(x){
  y <- 10*(x /(10+abs(x)))
  return(y)
}

# Make the plot
#bias.plot <- 
  ggplot()+
  geom_abline(slope = 1, intercept = 0, size = 1)+
  geom_abline(slope = 0.95, intercept = -0.5, linetype = 2, size = 1)+
  stat_function(aes(-20:20), fun = curve, n = 100, linetype = 3, size = 1)+
  lims(x = c(0,10), y = c(0,10))+
  labs(x = "True", y = "Estimated")+
  theme_classic(base_size = 18)

ggsave(bias.plot, file = "bias.jpeg", dpi = 650, width = 7, height = 5, units = "in")
