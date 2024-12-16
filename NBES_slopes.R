## packages ##
library(here)
library(tidyverse)
library(ggbeeswarm)

#### Load data ####
nbes_data100 <- readRDS("output/nbesSummary_press100.RData")
nbes_data100_flux <- readRDS("output/nbesSummary_fluctuation100.RData")
nbes_data100_combined <- readRDS("output/nbesSummary_combined100.RData")
communityMeta_press<-readRDS('output/nbesCommunityMeta_press100.RData')

#add diversity level description
nbes_plot_press <- nbes_data100 %>%
  mutate(RD = tOptUpper)

nbes_plot_press$RD[nbes_plot_press$RD==17.5] <- 'No RD (Topt = 17.5)'
nbes_plot_press$RD[nbes_plot_press$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_press$RD[nbes_plot_press$RD==20] <- 'High RD (15<Topt<20)'

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
  right_join(., nbes_plot_press)

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
ggsave(plot = last_plot(), file = here('output/exemplaryRun_press.png'))        


#### fluctuation ####
nbes_plot_flux <- nbes_data100_flux %>%
  mutate(RD = tOptUpper)

nbes_plot_flux$RD[nbes_plot_flux$RD==17.5] <- 'No RD (Topt = 15)'
nbes_plot_flux$RD[nbes_plot_flux$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_flux$RD[nbes_plot_flux$RD==20] <- 'High RD (15<Topt<20)'

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
  right_join(., nbes_plot_flux)
 
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
ggsave(plot = last_plot(), file = here('output/exemplaryRun_flux.png'))        


#### combined ####
#add diversity level description
nbes_plot_combined <- nbes_data100_combined %>%
  mutate(RD = tOptUpper)

nbes_plot_combined$RD[nbes_plot_combined$RD==17.5] <- 'No RD (Topt = 17.5)'
nbes_plot_combined$RD[nbes_plot_combined$RD==19] <- 'Med RD (16<Topt<19)'
nbes_plot_combined$RD[nbes_plot_combined$RD==20] <- 'High RD (15<Topt<20)'

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
  right_join(., nbes_plot_combined)

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
ggsave(plot = last_plot(), file = here('output/exemplaryRun_combined.png'))        


#### count ####
names(RD_press1)
negSlope_press <- RD_press1%>%  
  filter(resil.lm <0) %>%
  select(communityID, distType, tOptUpper, compNormSd, compNormMean,RD, resil.lm, combination, NBES, nSpecies)%>%
  group_by(communityID,nSpecies,RD,distType,compNormSd,compNormMean)%>%
  reframe(mean = mean(NBES))%>%
  spread(key = nSpecies, value = mean) %>%
  dplyr::rename(richness_2 = '2',
                richness_3 = '3',
                richness_4 = '4',
                richness_5 = '5')%>%
  filter(richness_5>richness_2) %>%
  pivot_longer(cols = c(richness_2, richness_3, richness_4, richness_5), names_to = 'richness_level', values_to  = 'mean.nbes')%>%
  separate(richness_level, into = c('removeme', 'nSpecies'))%>%
  mutate(nSpecies = as.numeric(nSpecies)) %>%
  left_join(., RD_press1)%>%
  select(nSpecies, combination, distType, communityID, NBES, mean.nbes,RD, resil.lm)%>%
  left_join(., communityMeta_press)

ggplot(negSlope_press, aes( x = nSpecies, y = NBES, color = RD))+
  geom_hline(yintercept = 0)+
  geom_smooth(se = F)+
  geom_point(size = 3)+
  facet_grid(~communityID)+
  theme_bw()+
  theme(legend.position = 'bottom')
ggsave(plot = last_plot(), file = here('output/V-shape.png'), width = 8, height = 4.5)


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

ggplot(GrandMean, aes( y = tOptValue, x = devFromGrandMean, color = nSpecies))  +
  geom_vline(xintercept = 0)+
  geom_point(size = 2)+
  labs(x = 'Influence on NBES', y='Topt', color = 'Treatment')+
  facet_wrap(~RD, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'bottom')
ggsave(plot=last_plot(), file = here('output/topt_NBESeffect-temp.png'), width = 8, height = 5)

  
#### other distType ####
negSlope_combined <- RD_combined1%>%  
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

unique(combined_negSlope$communityID)
press_negSlope <- RD_press1%>%  filter(resil.lm <0)
unique(press_negSlope$communityID)
fluct_negSlope <- RD_flux1%>%  filter(resil.lm <0)
unique(fluct_negSlope$communityID)
