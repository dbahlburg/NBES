simData <- readRDS('output/press_R5.RData')

simSettings <- simData[[1]]
simResults <- bind_rows(simData[[2]])

simMerged <- simResults %>% 
  left_join(., simSettings) %>% 
  filter(compNormSd %in% c(0, 0.5)) %>% 
  filter(tOptLower == 16 & tOptUpper == 19)

simIDs <- simMerged %>% 
  group_by(compNormSd) %>% 
  distinct(runID) %>% 
  mutate(runID = str_replace(runID, '_control',''),
         runID = str_replace(runID, '_treatment','')) %>% 
  group_by(compNormSd) %>% 
  distinct(runID) %>% 
  filter(str_detect(runID, 'S1S2S3S4S5')) %>% 
  slice(sample(1:20, 1))


simMerged %>% 
  filter(str_detect(runID, paste(c(simIDs$runID), collapse = '|'))) %>% 
  select(time, S1, S2, S3, S4, S5, compNormSd, distTreatment) %>% 
  pivot_longer(names_to = 'species', values_to = 'biomass', -c(time, compNormSd, distTreatment)) %>% 
  mutate(compNormSd = ifelse(compNormSd == 0, 'no competition', 'high competition')) %>% 
  ggplot(.,aes(x = time, y = biomass, colour = species)) +
  geom_line() +
  labs(x = 'time', y = 'biomass') +
  facet_grid(compNormSd~distTreatment) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = NA, colour = '#2b2b2b'), 
        panel.grid.major.y = element_line(colour = '#2b2b2b50', linetype = '22'),
        legend.position = 'bottom',
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 17),
        axis.title = element_text(size = 17),
        legend.title = element_text(size = 17))

ggsave('output/FigureS5.png', width = 10, height = 7)
ggsave('output/FigureS5.pdf', width = 10, height = 7)

  
  
  
  