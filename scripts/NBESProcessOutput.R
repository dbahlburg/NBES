# script to process an output file from the temperature-dependent community growth model.
# In the following script, we calculate the Net Biodiversity Effect on Stability (NBES)
# for communities with varying response diversities and degrees of competitiveness.
# ------------------------------------------------------------------------------------------ #
# load relevant packages
library(tidyverse)
library(MESS)
library(ggbeeswarm)
library(here)

# read model run for which NBES should be calculated
modelResults <- readRDS('output/press_R5.RData')
modelRuns <- modelResults[[2]]

# loop through model runs and extract runID
runIdentityTib <- NULL
for(k in 1:length(modelRuns)){
  runResults <- modelRuns[[k]]
  runIDTib <- tibble(itemPosition = k,
                     runID = unique(runResults$runID))
  runIdentityTib <- runIdentityTib %>% 
    bind_rows(., runIDTib)
}

# extract distinct communities that loop will iterate through
modelMeta <- modelResults[[1]]
distinctCommunities <- unique(modelMeta$communityID)

nbesDatAll <- NULL
communityMeta <- NULL
# iterate through communities and calculate NBES
for(i in 1:length(distinctCommunities)){
  
  # locate locations of list elements that correspond to the community
  listElementsInd <- runIdentityTib %>% 
    filter(str_detect(runID, paste(distinctCommunities[i],'_',sep='')))
  
  # we need to isolate the model runs with monocultures to calculate the expected Response Ratio:
  monoRunsMeta <- modelMeta %>% 
    filter(str_detect(runID, paste(distinctCommunities[i],'_',sep=''))) %>% 
    filter(nSpecies == 1)
  
  # first, we need to extract the biomass dynamics for each species in monoculture under treatment and control and merge everything in one table
  monocultures <- NULL
  speciesIDs <- unique(sapply(str_split(monoRunsMeta$runID, "_"), function(x) x[4]))
  for(j in 1:length(speciesIDs)){
    
    monoRunIDs <- listElementsInd %>% 
      filter(str_detect(runID, paste('_',speciesIDs[j],'_', sep='')))
    
    monoRuns <- bind_rows(modelRuns[c(monoRunIDs$itemPosition)]) %>% 
      pivot_wider(names_from = runID, values_from = speciesIDs[j]) %>% 
      rename(biomassMonoControl = 2, biomassMonoTreatment = 3) %>% 
      mutate(species = speciesIDs[j],
             biomassMonoRatio = biomassMonoTreatment/biomassMonoControl)
    
    monocultures <- monocultures %>% 
      bind_rows(., monoRuns)
  }
  
  # now we need to iterate through the model runs where the different species were simulated in mixed communities of varying richness
  # to make calculations more efficient, we calculate NBES separately for each richness level
  
  # get richness level:
  richnessLevels <- 2:max(modelMeta$maxSpecies)
  for(r in 1:length(richnessLevels)){
    
    # get metadata for model runs that correspond to current richness level
    mixedRunsMeta <- modelMeta %>% 
      filter(str_detect(runID, paste(distinctCommunities[i],'_',sep=''))) %>% 
      filter(nSpecies == richnessLevels[r])
    
    # iterate through different species combinations within that richness level
    speciesCombos <- unique(mixedRunsMeta$speciesCombo)
    for (m in 1:length(speciesCombos)){
      
      # identify model runs that correspond to current species combination
      speciesID <- unlist(speciesCombos[m])
      identical_elements <- which(sapply(mixedRunsMeta$speciesCombo, function(x) identical(x, speciesID)))
      currentRunMeta <- mixedRunsMeta %>% 
        slice(identical_elements)
      
      # now identify the corresponding elements in the list of model runs
      mixedRunIDs <- listElementsInd %>% 
        filter(runID %in% currentRunMeta$runID)
      
      # extract treatment and control runs for current species combination
      # create tidy tibble with control and treatment biomass for each species in mixed community
      mixedRun <- bind_rows(modelRuns[c(mixedRunIDs$itemPosition)]) %>% 
        pivot_longer(names_to = 'species', values_to = 'biomass', -c(time, runID)) %>% 
        pivot_wider(names_from = runID, values_from = biomass) %>% 
        rename(biomassMixControl = 3, biomassMixTreatment = 4) %>% 
        arrange(species, time)
      
      # create a second, tidy tibble, with total biomass in control and treatment in mixed community
      totalBiomassMixed <- mixedRun %>% 
        group_by(time) %>% 
        summarize(totalBiomMixControl = sum(biomassMixControl, na.rm = T),
                  totalBiomMixTreatment = sum(biomassMixTreatment, na.rm = T))
      
      # create object containing all variable required to calculate NBES
      masterDat <- mixedRun %>% 
        left_join(., monocultures, by = join_by(time, species)) %>% 
        left_join(., totalBiomassMixed, by = join_by(time)) %>% 
        ungroup()%>%
        mutate(combination = paste('S',speciesID,collapse = '', sep = ''),
               relBiomassT0Mixed = unique(mixedRunsMeta$N0)/(unique(mixedRunsMeta$N0) * unique(mixedRunsMeta$nSpecies))) %>% 
        select(time, combination, species, biomassMixTreatment, biomassMixControl,
               biomassMonoTreatment,biomassMonoControl, biomassMonoRatio, totalBiomMixTreatment, totalBiomMixControl,relBiomassT0Mixed) %>%
        mutate(expSppBiom = totalBiomMixControl*biomassMonoRatio*relBiomassT0Mixed,
               RRexp = ((expSppBiom-relBiomassT0Mixed*totalBiomMixControl)/
                        (expSppBiom+relBiomassT0Mixed*totalBiomMixControl))) %>%
        group_by(time, combination) %>%
        mutate(RRobs = (totalBiomMixTreatment-totalBiomMixControl)/(totalBiomMixTreatment+totalBiomMixControl),
               SumExpSppBiom = sum(expSppBiom, na.rm = T),
               RR_ges_exp = (SumExpSppBiom-totalBiomMixControl)/(SumExpSppBiom+totalBiomMixControl))
      
      
      # extract competition matrix to summarise competition metrics about community
      compMatrixMixedRun <- do.call(rbind, currentRunMeta$compMatrix[1])[speciesID,speciesID]

      # calculate nbes, add data about competitiveness of community
      nbesDat <- masterDat %>%
        group_by(combination) %>%
        summarise(AUC.RR_obs= auc(time, RRobs,  from = min(time, na.rm = TRUE), to = max(150, na.rm = TRUE),
                                  type = c("linear"),absolutearea = FALSE),
                  AUC.RR_spp_exp= auc(time, RRexp,  from = min(time, na.rm = TRUE), to = max(150, na.rm = TRUE),
                                  type = c("linear"),absolutearea = FALSE),
                  AUC.RR_exp=auc(time, RR_ges_exp,  from = min(time, na.rm = TRUE), to = max(150, na.rm = TRUE),
                                 type = c("linear"),absolutearea = FALSE),
                  NBES = AUC.RR_obs-AUC.RR_exp) %>% 
        mutate(speciesCombo = list(speciesID),
               meanAlphas = mean(c(compMatrixMixedRun[upper.tri(compMatrixMixedRun)], compMatrixMixedRun[lower.tri(compMatrixMixedRun)])),
               sdAlphas = sd(c(compMatrixMixedRun[upper.tri(compMatrixMixedRun)], compMatrixMixedRun[lower.tri(compMatrixMixedRun)]))) %>% 
        left_join(., filter(mixedRunsMeta, str_detect(runID, 'treatment')), by = 'speciesCombo') %>% 
        select(-compMatrix)
      
      nbesDatAll <- nbesDatAll %>% 
        bind_rows(., nbesDat)
      
      # store species-specific alpha values of current model run/community
      for (f in 1:length(speciesID)){
        alphaValues <- nbesDat %>% 
          select(runID, communityID, combination, tOptValues) %>% 
          mutate(species = paste('spp', speciesID[f], sep = ''),
                 meanAlpha = mean(compMatrixMixedRun[,f][-f]),
                 tOptValue = unlist(tOptValues)[f]
                 ) %>% 
          select(-tOptValues)
        communityMeta <- communityMeta %>% 
          bind_rows(., alphaValues)
      }
      
      
      }
    
  }
  print(paste('community', i, 'out of', length(distinctCommunities),'processed', sep = ' '))
}

# save summary file
write_rds(nbesDatAll, 'output/nbesSummary_press150.RData')
write_rds(communityMeta, 'output/nbesCommunityMeta_press150.RData')

nbesDatAll %>% 
  ggplot(.,aes(x = nSpecies, y = NBES)) +
  geom_hline(yintercept = 0)+
  geom_quasirandom()


### one exemplary run ###
write_rds(masterDat, 'output/masterDat_press.RData')

 masterDat %>% 
   ggplot(.,aes(x = time, y = totalBiomMixControl)) +
   geom_line(color ='black')+
   geom_line(aes(y = totalBiomMixTreatment), color = 'darkred')+
   labs(y = 'Total Biomass')
 ggsave(plot = last_plot(), file = here('output/TotalBiomass_pressDist.png'))

 masterDat %>% 
   ggplot(.,aes(x = time, y = biomassMonoTreatment, color = species)) +
   geom_line()+
   geom_line(aes(y = biomassMonoControl), linetype = 'dashed')+
   labs(y = 'Total Biomass')
 ggsave(plot = last_plot(), file = here('output/sppBiomass_combinedDist.png'))
 




 
 






