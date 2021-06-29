# Creating example environmental covariate figure
# Part of fig. S1; prior augmentation manuscript
# Summer 2021

# Load package
library(tidyverse)

# Create plot
ggplot()+
  geom_abline(slope = 0.6, intercept = 1, size = 1)+
  scale_x_continuous(limits = c(0,5), expand = c(0,0))+
  scale_y_continuous(limits = c(0,5), expand = c(0,0))+
  labs(x = expression(beta[1]), y = expression(psi[ij]))+
  theme_bw(base_size = 42)+
  theme(panel.grid = element_blank(), axis.text = element_blank())

ggsave(filename = "samplecov.jpeg", dpi = 600)
