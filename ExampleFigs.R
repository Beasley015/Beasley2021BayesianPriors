# Creating example figures
# Part of fig. S1; prior augmentation manuscript
# Summer 2021

# Load package
library(tidyverse)

# Sample occupancy distribution
ggplot()+
  stat_function(fun = dbeta, n = 100, args = list(shape1 = 7, shape2 = 4),
                size = 1)+
  # stat_function(fun = dbeta, n = 100, args = list(shape1 = 7, shape2 = 15), 
  #               color = "red", size = 1)+
  # stat_function(fun = dbeta, n = 100, args = list(shape1=10, shape2=10),
  #               color = "red", linetype = "dashed", size = 1)+
  scale_x_continuous(breaks = c(0, 1))+
  labs(x = expression(psi))+
  theme_bw(base_size = 36)+
  theme(panel.grid = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(margin = margin(t = -30, r = 0, b = 0, 
                                                    l = 0)))

ggsave(filename = 'sampledist.jpeg', dpi = 600, width = 6, height = 4)

# Create sample environmental covariate figure
ggplot()+
  geom_abline(slope = 0.6, intercept = 1, size = 1)+
  scale_x_continuous(limits = c(0,5), expand = c(0,0))+
  scale_y_continuous(limits = c(0,5), expand = c(0,0))+
  labs(x = expression(beta[1]), y = expression(logit(psi[ij])))+
  theme_bw(base_size = 42)+
  theme(panel.grid = element_blank(), axis.text = element_blank())

ggsave(filename = "samplecov.jpeg", dpi = 600, width = 6, height = 4)
