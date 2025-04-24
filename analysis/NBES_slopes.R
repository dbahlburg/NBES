## packages ##
library(here)
library(tidyverse)
library(ggbeeswarm)

#### Load data ####
nbes_data100 <- readRDS("output/nbesSummary_press100.RData")
nbes_data100_flux <- readRDS("output/nbesSummary_fluctuation100.RData")
nbes_data100_combined <- readRDS("output/nbesSummary_combined100.RData")
communityMeta_press<-readRDS('output/nbesCommunityMeta_press100.RData')
communityMeta_combined<-readRDS('output/nbesCommunityMeta_combined100.RData')
communityMeta_fluctuation<-readRDS('output/nbesCommunityMeta_fluctuation100.RData')

#### press ####
#add diversity level description
nbes_plot_press <- nbes_data100 %>%
  mutate(RD = tOptUpper)

nbes_plot_press$RD[nbes_plot_press$RD==17.5] <- 'same TOpt'
nbes_plot_press$RD[nbes_plot_press$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_press$RD[nbes_plot_press$RD==20] <- 'different TOpt'

names(nbes_plot_press)
USI_press <- unique(nbes_plot_press$communityID)

#create an empty data frame
slope_press<-tibble()

# the following loop cycles through all unique cases

for(i in 1:length(USI_press)){
  temp<-nbes_plot_press[nbes_plot_press$communityID==USI_press[i], ]#creates a temporary data frame for each case
  lm1<-lm(NBES~log(nSpecies), temp)#makes a linear regression
  intcp.lm <- coef(summary(lm1))[1, 1]#selects the intercept
  se.intcp.lm<- coef(summary(lm1))[1, 2]#selects its standard error
  resil.lm <- coef(summary(lm1))[2, 1]#selects the slope
  se.slp.lm<- coef(summary(lm1))[2, 2]#selects its standard error
  sd.res.lm<- sd(resid(lm1)) #selects the standard deviation of the residuals
  slope_press<-rbind(slope_press,data.frame(temp[1,c(10, 11,14, 17, 18,32)],resil.lm,intcp.lm,se.intcp.lm,sd.res.lm,se.slp.lm))
  rm(temp)
}
names(slope_press)


hist(slope_press$resil.lm)
ggplot(slope_press, aes(x = compNormSd, y = resil.lm))+
  geom_point()+
  facet_grid(~RD)

RD_press1 <- slope_press%>%
  right_join(., nbes_plot_press)%>%
  mutate(alpha_info = paste('alpha','==', compNormSd))

RD_press1%>%
  filter(resil.lm <0) %>%
  group_by(nSpecies, RD, compNormSd, communityID) %>%
  reframe(mean = mean(NBES))%>%
  ggplot(., aes (x = nSpecies, y = mean, color = communityID))+
  geom_hline(yintercept = 0)+
  geom_quasirandom()+
  geom_line()+
  facet_grid(~RD~compNormSd, scales = 'free_y')+
  theme(legend.position = 'none')
#ggsave(plot = last_plot(), file = here('output/exemplaryRun_press.png'))        


#### fluctuation ####
nbes_plot_flux <- nbes_data100_flux %>%
  mutate(RD = tOptUpper)

nbes_plot_flux$RD[nbes_plot_flux$RD==17.5] <- 'same TOpt'
nbes_plot_flux$RD[nbes_plot_flux$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_flux$RD[nbes_plot_flux$RD==20] <- 'different TOpt'

names(nbes_plot_flux)
USI_flux <- unique(nbes_plot_flux$communityID)

#create an empty data frame
slope_fluct<-tibble()

# the following loop cycles through all unique cases

for(i in 1:length(USI_flux)){
  temp<-nbes_plot_flux[nbes_plot_flux$communityID==USI_flux[i], ]#creates a temporary data frame for each case
  lm1<-lm(NBES~log(nSpecies), temp)#makes a linear regression
   intcp.lm <- coef(summary(lm1))[1, 1]#selects the intercept
   se.intcp.lm<- coef(summary(lm1))[1, 2]#selects its standard error
   resil.lm <- coef(summary(lm1))[2, 1]#selects the slope
   se.slp.lm<- coef(summary(lm1))[2, 2]#selects its standard error
   sd.res.lm<- sd(resid(lm1)) #selects the standard deviation of the residuals
   slope_fluct<-rbind(slope_fluct,data.frame(temp[1,c(10, 11,14, 17, 18,32)],resil.lm,intcp.lm,se.intcp.lm,sd.res.lm,se.slp.lm))
  rm(temp)
}
names(slope_fluct)


hist(slope_fluct$resil.lm)
ggplot(slope_fluct, aes(x = compNormSd, y = resil.lm))+
  geom_point()+
  facet_grid(~RD)

RD_flux1 <- slope_fluct%>%
  right_join(., nbes_plot_flux)%>%
  mutate(alpha_info = paste('alpha','==', compNormSd))
 
RD_flux1%>%
  filter(resil.lm <0) %>%
  group_by(nSpecies, RD, compNormSd, communityID) %>%
  reframe(mean = mean(NBES))%>%
  ggplot(., aes (x = nSpecies, y = mean, color = communityID))+
  geom_hline(yintercept = 0)+
  geom_quasirandom()+
  geom_line()+
  facet_grid(~RD~compNormSd, scales = 'free_y')+
  theme(legend.position = 'none')
#ggsave(plot = last_plot(), file = here('output/exemplaryRun_flux.png'))        


#### combined ####
#add diversity level description
nbes_plot_combined <- nbes_data100_combined %>%
  mutate(RD = tOptUpper)

nbes_plot_combined$RD[nbes_plot_combined$RD==17.5] <- 'same TOpt'
nbes_plot_combined$RD[nbes_plot_combined$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_combined$RD[nbes_plot_combined$RD==20] <- 'different TOpt'

names(nbes_plot_combined)
USI_combined <- unique(nbes_plot_combined$communityID)

#create an empty data frame
slope_combined<-tibble()

# the following loop cycles through all unique cases

for(i in 1:length(USI_combined)){
  temp<-nbes_plot_combined[nbes_plot_combined$communityID==USI_combined[i], ]#creates a temporary data frame for each case
  lm1<-lm(NBES~log(nSpecies), temp)#makes a linear regression
  intcp.lm <- coef(summary(lm1))[1, 1]#selects the intercept
  se.intcp.lm<- coef(summary(lm1))[1, 2]#selects its standard error
  resil.lm <- coef(summary(lm1))[2, 1]#selects the slope
  se.slp.lm<- coef(summary(lm1))[2, 2]#selects its standard error
  sd.res.lm<- sd(resid(lm1)) #selects the standard deviation of the residuals
  slope_combined<-rbind(slope_combined,data.frame(temp[1,c(10, 11,17,18,32)],resil.lm,intcp.lm,se.intcp.lm,sd.res.lm,se.slp.lm))
  rm(temp)
}
names(slope_combined)


hist(slope_combined$resil.lm)
ggplot(slope_combined, aes(x = compNormSd, y = resil.lm))+
  geom_point()+
  facet_grid(~RD)

RD_combined1 <- slope_combined%>%
  right_join(., nbes_plot_combined)%>%
  mutate(alpha_info = paste('alpha','==', compNormSd))

RD_combined1%>%
  filter(resil.lm <0) %>%
  group_by(nSpecies, RD, compNormSd, communityID) %>%
  reframe(mean = mean(NBES))%>%
  ggplot(., aes (x = nSpecies, y = mean, color = communityID))+
  geom_hline(yintercept = 0)+
  geom_quasirandom()+
  geom_line()+
  facet_grid(~RD~compNormSd, scales = 'free_y')+
  theme(legend.position = 'none')
#ggsave(plot = last_plot(), file = here('output/exemplaryRun_combined.png'))        

#### plot ####
str(RD_press1 )

press_plot <- RD_press1%>%
  group_by(nSpecies, RD, alpha_info, communityID) %>%
  filter( tOptUpper == 20)%>%
  reframe(mean = mean(NBES),
          sd = sd(NBES),
          se = sd/sqrt(n()))%>%
  ggplot(., aes (x = nSpecies, y = mean, group = communityID))+
  geom_hline(yintercept = 0)+
 # geom_errorbar(aes(ymin = mean-se, ymax = mean+se))+ 
  geom_point(alpha = 0.5, color = "#0072B2")+
  geom_line(alpha = 0.6, color = "#0072B2")+
  labs(y = 'Mean NBES', title = 'Tempterature Increase', x = 'Species Richness')+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
press_plot

flux_plot <- RD_flux1%>%
  group_by(nSpecies, RD, alpha_info, communityID) %>%
  filter( tOptUpper == 20)%>%
  reframe(mean = mean(NBES),
          sd = sd(NBES),
          se = sd/sqrt(n()))%>%
  ggplot(., aes (x = nSpecies, y = mean, group = communityID))+
  geom_hline(yintercept = 0)+
  # geom_errorbar(aes(ymin = mean-se, ymax = mean+se))+ 
  geom_point(alpha = 0.5, color = "#E41A1C")+
  geom_line(alpha = 0.6, color = "#E41A1C")+
  labs(y = 'Mean NBES', title = 'Temperature Fluctuations', x = 'Species Richness')+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))

flux_plot

combined_plot <- RD_combined1%>%
  group_by(nSpecies, RD, alpha_info, communityID) %>%
  filter( tOptUpper ==20)%>%
  reframe(mean = mean(NBES),
          sd = sd(NBES),
          se = sd/sqrt(n()))%>%
  ggplot(., aes (x = nSpecies, y = mean, group = communityID))+
  geom_hline(yintercept = 0)+
  # geom_errorbar(aes(ymin = mean-se, ymax = mean+se))+ 
  geom_point(alpha = 0.5, color = '#c7b514')+
  geom_line(alpha = 0.6, color = '#c7b514')+
  labs(y = 'Mean NBES', title = 'Temperature Increase + Fluctuations', x = 'Species Richness')+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
combined_plot

cowplot::plot_grid(press_plot, flux_plot, combined_plot, labels = c('(a)', '(b)', '(c)'), ncol = 1)
ggsave(plot = last_plot(), file = here('output/ExtendedData_FigureS1_MeanNBES_allDist_RD.png'), width = 8, height = 9)

press_plot1 <- RD_press1%>%
  group_by(nSpecies, RD, alpha_info, communityID) %>%
  filter( tOptUpper ==17.5)%>%
  reframe(mean = mean(NBES),
          sd = sd(NBES),
          se = sd/sqrt(n()))%>%
  ggplot(., aes (x = nSpecies, y = mean, group = communityID))+
  geom_hline(yintercept = 0)+
  # geom_errorbar(aes(ymin = mean-se, ymax = mean+se))+ 
  geom_point(alpha = 0.6, color = "#0072B2")+
  geom_line(alpha = 0.5, color = "#0072B2")+  
  labs(y = 'Mean NBES', title = 'Same Temperature Optima', x = 'Species Richness')+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
press_plot1

press_plot2 <- RD_press1%>%
  group_by(nSpecies, RD, alpha_info, communityID) %>%
  filter( tOptUpper ==20)%>%
  reframe(mean = mean(NBES),
          sd = sd(NBES),
          se = sd/sqrt(n()))%>%
  ggplot(., aes (x = nSpecies, y = mean, group = communityID))+
  geom_hline(yintercept = 0)+
  # geom_errorbar(aes(ymin = mean-se, ymax = mean+se))+ 
  geom_point(alpha = 0.6, color = "#0072B2")+
  geom_line(alpha = 0.5, color = "#0072B2")+
  labs(y = 'Mean NBES', title = 'Different Tempterature Optima', x = 'Species Richness')+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 10,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 10,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.background =element_rect(),
        strip.text.x  = element_text(size = 10, face = 'plain'))+
  guides(color = guide_legend(override.aes = list(size = 3.5)))+
  theme(legend.position = 'none',
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14))
press_plot2

cowplot::plot_grid(press_plot1, press_plot2, labels = c('(a)', '(b)', '(c)'), ncol = 1)
ggsave(plot = last_plot(), file = here('output/Figure2_MeanNBES_press.pdf'), width = 7, height = 5)

#### Grand Mean ####
names(RD_press1)
negSlope_press <- RD_press1%>%  
  #filter(resil.lm <0) %>%
  select(communityID, distType, tOptUpper, compNormSd, alpha_info,compNormMean,RD, resil.lm, combination, NBES, nSpecies)%>%
  group_by(communityID,nSpecies,RD,distType,compNormSd,alpha_info,compNormMean)%>%
  reframe(mean = mean(NBES))%>%
  spread(key = nSpecies, value = mean) %>%
  dplyr::rename(richness_2 = '2',
                richness_3 = '3',
                richness_4 = '4',
                richness_5 = '5')%>%
 # filter(richness_2>richness_4) %>%
  pivot_longer(cols = c(richness_2, richness_3, richness_4, richness_5), names_to = 'richness_level', values_to  = 'mean.nbes')%>%
  separate(richness_level, into = c('removeme', 'nSpecies'))%>%
  mutate(nSpecies = as.numeric(nSpecies)) %>%
  left_join(., RD_press1)%>%
  select(nSpecies, combination, distType, communityID, NBES, mean.nbes,RD, alpha_info,compNormSd,resil.lm)%>%
  left_join(., communityMeta_press)

ggplot(negSlope_press, aes( x = nSpecies, y = NBES, color = RD))+
  geom_hline(yintercept = 0)+
  geom_point(size = 3)+
  facet_grid(~alpha_info, scales = 'free_y', labeller = label_parsed)+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot = last_plot(), file = here('output/V-shape.png'), width = 8, height = 4.5)


S1_negSlope_press <- negSlope_press %>%
  filter(str_detect(combination, 'S1')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp1'))

S2_negSlope_press <- negSlope_press %>%
  filter(str_detect(combination, 'S2')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp2'))

S3_negSlope_press <- negSlope_press %>%
  filter(str_detect(combination, 'S3')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp3'))

S4_negSlope_press <- negSlope_press %>%
  filter(str_detect(combination, 'S4')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp4'))

S5_negSlope_press <- negSlope_press %>%
  filter(str_detect(combination, 'S5')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp5'))

allS<- S5_negSlope_press %>%
  bind_rows(., S4_negSlope_press) %>%
  bind_rows(., S3_negSlope_press) %>%
  bind_rows(., S2_negSlope_press) %>%
  bind_rows(., S1_negSlope_press) 

GrandMean <- negSlope_press %>%
  group_by(nSpecies, communityID,RD)%>%
  reframe(gMean = mean(NBES)) %>%
  right_join(., allS) %>%
  mutate(devFromGrandMean = S_nbes-gMean) %>%
  left_join(., negSlope_press)


#set colour palette 
colours<-c("darkblue", "skyblue", '#EFC000FF', 'darkorange')
custom_labeller <- labeller(
  alpha_info = label_parsed,    # this one gets parsed
  group = label_value      # this one stays plain
)
ggplot(GrandMean, aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha = 0.4)+
  labs(x = 'Influence on NBES',  color = 'Richness')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~RD~alpha_info,scales = 'free_y', labeller=custom_labeller)+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot=last_plot(), file = here('output/topt_NBESeffect-temp.png'), width = 8, height = 9)

  
#### other distType ####
negSlope_combined <- RD_combined1%>%  
  filter(resil.lm <0) %>%
  select(communityID, distType, tOptUpper, compNormSd, alpha_info,RD, resil.lm, combination, NBES, nSpecies)%>%
  group_by(communityID,nSpecies,RD,distType,alpha_info,compNormSd)%>%
  reframe(mean = mean(NBES))%>%
  spread(key = nSpecies, value = mean) %>%
  dplyr::rename(richness_2 = '2',
                richness_3 = '3',
                richness_4 = '4',
                richness_5 = '5')%>%
 # filter(richness_2>richness_4)%>%
  pivot_longer(cols = c(richness_2, richness_3, richness_4, richness_5), names_to = 'richness_level', values_to  = 'mean.nbes')%>%
  separate(richness_level, into = c('removeme', 'nSpecies'))%>%
  mutate(nSpecies = as.numeric(nSpecies)) %>%
  left_join(., RD_combined1)%>%
  select(nSpecies, combination, distType, communityID, NBES, alpha_info,compNormSd,mean.nbes,RD, resil.lm)%>%
  left_join(., communityMeta_combined)

S1_negSlope_combined <- negSlope_combined %>%
  filter(str_detect(combination, 'S1')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp1'))

S2_negSlope_combined <- negSlope_combined %>%
  filter(str_detect(combination, 'S2')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp2'))

S3_negSlope_combined <- negSlope_combined %>%
  filter(str_detect(combination, 'S3')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp3'))

S4_negSlope_combined <- negSlope_combined %>%
  filter(str_detect(combination, 'S4')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp4'))

S5_negSlope_combined <- negSlope_combined %>%
  filter(str_detect(combination, 'S5')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp5'))

allS_combined<- S5_negSlope_combined %>%
  bind_rows(., S4_negSlope_combined) %>%
  bind_rows(., S3_negSlope_combined) %>%
  bind_rows(., S2_negSlope_combined) %>%
  bind_rows(., S1_negSlope_combined) 

GrandMean_combined <- negSlope_combined %>%
  group_by(nSpecies, communityID,RD)%>%
  reframe(gMean = mean(NBES)) %>%
  right_join(., allS_combined) %>%
  mutate(devFromGrandMean = S_nbes-gMean) %>%
  left_join(., negSlope_combined)

GrandMean_combined%>%
  filter(RD!="same TOpt" )%>%
ggplot(., aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES',  color = 'Richness')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~RD~alpha_info, labeller = custom_labeller, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot=last_plot(), file = here('output/topt_NBESeffect-temp_combined.png'), width = 8, height = 8)


###
negSlope_flux <- RD_flux1%>%  
  filter(resil.lm <0) %>%
  select(communityID, distType, tOptUpper, alpha_info, compNormSd, RD, resil.lm, combination, NBES, nSpecies)%>%
  group_by(communityID,nSpecies,RD,distType,alpha_info, compNormSd)%>%
  reframe(mean = mean(NBES))%>%
  spread(key = nSpecies, value = mean) %>%
  dplyr::rename(richness_2 = '2',
                richness_3 = '3',
                richness_4 = '4',
                richness_5 = '5')%>%
 # filter(richness_2>richness_4)%>%
  pivot_longer(cols = c(richness_2, richness_3, richness_4, richness_5), names_to = 'richness_level', values_to  = 'mean.nbes')%>%
  separate(richness_level, into = c('removeme', 'nSpecies'))%>%
  mutate(nSpecies = as.numeric(nSpecies)) %>%
  left_join(., RD_flux1)%>%
  select(nSpecies, combination, distType, communityID, NBES, alpha_info,compNormSd,mean.nbes,RD, resil.lm)%>%
  left_join(., communityMeta_fluctuation)

S1_negSlope_flux <- negSlope_flux %>%
  filter(str_detect(combination, 'S1')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp1'))

S2_negSlope_flux <- negSlope_flux %>%
  filter(str_detect(combination, 'S2')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp2'))

S3_negSlope_flux <- negSlope_flux %>%
  filter(str_detect(combination, 'S3')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp3'))

S4_negSlope_flux <- negSlope_flux %>%
  filter(str_detect(combination, 'S4')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp4'))

S5_negSlope_flux <- negSlope_flux %>%
  filter(str_detect(combination, 'S5')) %>%
  group_by(communityID, nSpecies) %>%
  reframe(S_nbes = mean(NBES),
          species = paste('spp5'))

allS_flux<- S5_negSlope_flux %>%
  bind_rows(., S4_negSlope_flux) %>%
  bind_rows(., S3_negSlope_flux) %>%
  bind_rows(., S2_negSlope_flux) %>%
  bind_rows(., S1_negSlope_flux) 

GrandMean_flux <- negSlope_flux %>%
  group_by(nSpecies, communityID,RD)%>%
  reframe(gMean = mean(NBES)) %>%
  right_join(., allS_flux) %>%
  mutate(devFromGrandMean = S_nbes-gMean) %>%
  left_join(., negSlope_flux)

GrandMean_flux%>%
  #filter(RD=="No RD (Topt = 17.5)" )%>%
  ggplot(., aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES',  color = 'Richness')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~RD~alpha_info, labeller = custom_labeller, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')
#ggsave(plot=last_plot(), file = here('output/topt_NBESeffect-temp_flux.png'), width = 8, height = 8)

#### Topt Plot ####
label_alpha <- c('no interaction', 'intermed interaction' , 'high interaction' ) #compNormSd
pa <- GrandMean %>%
  filter(RD=="different TOpt" )%>%
  ggplot(., aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha = 0.4)+
  labs(x = 'Influence on NBES',  color = 'Richness', title = 'Increase')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~alpha_info, labeller = label_parsed,scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))

ggsave(plot = pa, file = here('output/GrandMean_increase.png'), width = 8, height = 5)

pb<-GrandMean_combined %>%
  filter(RD=="different TOpt" )%>%
  ggplot(., aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha = 0.4)+
  labs(x = 'Influence on NBES', color = 'Richness', title = 'Increase + Fluctuation')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~alpha_info, labeller = label_parsed,scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))


pc<-GrandMean_flux %>%
  filter(RD=="different TOpt" )%>%
  ggplot(., aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha = 0.4)+
  labs(x = 'Influence on NBES', color = 'Richness', title = 'Fluctuation')+
  scale_color_gradientn(colours = colours)+
  ylab(bquote(b[opt]))+
  facet_grid(~alpha_info, labeller = label_parsed,scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')+
  theme(axis.title.x = element_text(size = 14,face = "plain", colour = "black", vjust = 0),
        axis.text.x = element_text(size = 12,  colour = "black", angle = 0, vjust = 0.5)) +
  theme(axis.title.y = element_text(size = 14, face = "plain", colour = "black", vjust = 1.8),
        axis.text.y = element_text(size = 12,  colour = "black", angle = 0, hjust = 0.4)) +
  theme(strip.text.x  = element_text(size = 12))


legend <- cowplot::get_legend(pa)
cowplot::plot_grid(pa+theme(legend.position = 'none'),
                   pc+theme(legend.position = 'none'),
                   pb+theme(legend.position = 'none'),
                   labels = c('(a)', '(b)','(c)'),ncol =1)
ggsave(plot=last_plot(), file = here('output/Figure3_Topt_GrandMean.png'), width = 7, height = 8)


#### counts #####
negSlope_fluct <- RD_flux1%>%  
  filter(resil.lm <0) %>%
  select(communityID, distType, tOptUpper, compNormSd, RD, resil.lm, combination, NBES, nSpecies)%>%
  group_by(communityID,nSpecies,RD,distType,compNormSd)%>%
  reframe(mean = mean(NBES))%>%
  spread(key = nSpecies, value = mean) %>%
  dplyr::rename(richness_2 = '2',
                richness_3 = '3',
                richness_4 = '4',
                richness_5 = '5')%>%
  filter(richness_5>richness_2)

