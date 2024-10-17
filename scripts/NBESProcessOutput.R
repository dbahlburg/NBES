# script to process an output file from the temperature-dependent community growth model.
# In the following script, we calculate the Net Biodiversity Effect on Stability (NBES)
# for communities with varying response diversities and degrees of competitiveness.
# ------------------------------------------------------------------------------------------ #
# load relevant packages
library(tidyverse)

# read model run for which NBES should be calculated
modelResults <- readRDS('/Users/dobahl001/github/NBES/output/press_R5.RData')
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
      rename(biomassControl = 2, biomassTreatment = 3) %>% 
      mutate(species = speciesIDs[j],
             biomassRatio = biomassTreatment/biomassControl)
    
    monocultures <- monocultures %>% 
      bind_rows(., monoRuns)
  }
  
  # now we need to iterate through the model runs where the different species were simulated in mixed communities of varying richness
  mixedRunsMeta <- modelMeta %>% 
    filter(str_detect(runID, paste(distinctCommunities[i],'_',sep=''))) %>% 
    filter(nSpecies > 1)
  
  # iterate through different species combinations
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
      rename(biomassControl = 3, biomassTreatment = 4) %>% 
      arrange(species, time)
    
    # create a second, tidy tibble, with total biomass in control and treatment in mixed community
      totalBiomassMixed <- mixedRun %>% 
        group_by(time) %>% 
        summarize(totalBiomControl = sum(biomassControl, na.rm = T),
                  totalBiomTreatment = sum(biomassTreatment, na.rm = T))
      
  }
  
}





















