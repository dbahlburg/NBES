### R script to analyse output ###
library(here)
library(tidyverse)

nbes_data100 <- readRDS('output/nbesSummary_100.RData')%>%
  mutate(RD = paste('high'))
#nbes_data <- readRDS('output/nbesSummary_highRD_sd0.RData')%>%
#  mutate(RD = paste('high'))

#add diversity level description
nbes_plot <- nbes_data100 %>%
  mutate(RD = tOptUpper)

nbes_plot$RD[nbes_plot$RD==15] <- 'No RD (Topt = 15)'
nbes_plot$RD[nbes_plot$RD==18] <- 'Low RD (17<Topt<18)'
nbes_plot$RD[nbes_plot$RD==19] <- 'Intermed RD (16<Topt<19)'
nbes_plot$RD[nbes_plot$RD==20] <- 'High RD (15<Topt<20)'

str(nbes_data)

### create plots ###

# NBES - Richness #
nbes_plot %>%
  group_by(nSpecies, RD, compNormSd)%>%
  mutate(mean.NBES = mean(NBES),
         sd.NBES = sd(NBES))%>%
  ggplot(., aes(x= nSpecies, y = NBES))+
    geom_hline(yintercept = 0)+
    geom_quasirandom(size = 0.7, alpha = 0.3)+
    geom_errorbar(aes(ymin = mean.NBES-sd.NBES, ymax = mean.NBES+sd.NBES), width = .1, color = 'black')+
    geom_point(aes(y = mean.NBES), color = 'darkred')+
   scale_x_continuous(limits = c(1.5,5.5),breaks = seq(2,5,1))+
    facet_grid(RD~compNormSd)+
    theme_bw()+
    theme(legend.position = 'none')
ggsave(plot = last_plot(), width = 8, height = 8, file = here('output/NBES_richness.png'))

# Competitive communities - NBES #
nbes_plot %>%
  filter(compNormSd ==0.5)%>%
  ggplot(., aes(x= meanAlphas, y = NBES))+
  geom_hline(yintercept = 0)+
  labs(x = 'Mean Community Alpha')+
  facet_grid(RD~nSpecies)+
  geom_point(size = 0.7, alpha = 0.4)+
  geom_smooth(method = 'lm')+
  theme_bw()+
  theme(legend.position = 'none')
ggsave(plot = last_plot(), width = 8, height = 8, file = here('output/NBES_MeanAlpha.png'))

# No competition - NBES #
nbes_plot %>%
  filter(compNormSd == 0)%>%
  ggplot(., aes(x= nSpecies, y = NBES))+
  facet_grid(~RD)+
  geom_hline(yintercept = 0)+
  geom_quasirandom(size = 0.7, alpha = 0.3)+
  theme_bw()
ggsave(plot = last_plot(), width = 8, height = 4, file = here('output/NBES_Richness_sd0.png'))

