runScenarios <- function(destPath, maxSpecies, distType, compType, pTempMin = NA, pTempMax = NA, tempControl, fTempMin = NA, fTempMax = NA, fTempPeriod = NA, compNormMean, repetitions, N0, beta, delta, a_b, s, a_d, z, tmax, dt, compNormSd, tOptTib, saveFile = T){
  
  # load auxiliary functions to run the model
  source('scripts/NBESGenerateCompetitionMatrix.R')
  source('scripts/NBESGenerateTOpt.R')
  source('scripts/NBESGenerateTemperatures.R')
  source('scripts/NBESrk4Solver.R')
  source('scripts/NBESGrowthRates.R')
  
  # calculate possible species combinations and store in list:
  speciesCombTib <- NULL
  for(i in 1:maxSpecies){
    speciesComb <- t(combn(1:maxSpecies, i))
    for(rowInd in 1:nrow(speciesComb)){
      
      speciesCombSubTib <- tibble(nSpecies = length(speciesComb[rowInd,]),
                                  speciesCombo = list(speciesComb[rowInd,]))
      
      speciesCombTib <- speciesCombTib %>% 
        bind_rows(., speciesCombSubTib)
    }
  }
  
  # create competition matrices (for each compNormSd value + 10 repetitions to account for stochasticity)
  compMatrixTib <- NULL
  for(k in 1:length(compNormSd)){
    for(j in 1:repetitions){
      compMatrix <- generateCompetitionMatrix(competitionScenario = compType, 
                                              nSpecies = maxSpecies,
                                              meanNormComp = compNormMean,
                                              sdNormComp = compNormSd[k])
      compMatrixSubTib <- tibble(compNormSd = compNormSd[k],
                                 repetition = j,
                                 compMatrix = list(compMatrix))
      compMatrixTib <- compMatrixTib %>% 
        bind_rows(.,compMatrixSubTib)
    }
  }
  
  # create temperature time series
  tempTib <- tibble(distTreatment = c('control','treatment'),
                    temperatureValues = c(list(rep(tempControl, times = tmax/dt)),
                                          list(generateTemperatures(distType = distType, 
                                                                    pTempMin = pTempMin, 
                                                                    pTempMax = pTempMax, 
                                                                    fTempMin = fTempMin, 
                                                                    fTempMax = fTempMax, 
                                                                    fTempPeriod = fTempPeriod, 
                                                                    tmax = tmax, dt = dt))))
                      
                      
  # create tibble with all simulation settings
  simulationSettingsTib <- tOptTib %>% 
    crossing(compMatrixTib) %>%
    crossing(speciesCombTib) %>%
    crossing(tempTib) %>% 
    mutate(runID = row_number(),
           distType = distType,
           maxSpecies = maxSpecies,
           compType = compType,
           compNormMean = compNormMean, 
           N0 = N0,
           beta = beta,
           delta = delta,
           a_b = a_b,
           s = s,
           a_d = a_d,
           z = z,
           tmax = tmax,
           dt = dt) %>% 
    select(runID, distType, distTreatment, compType, compNormMean, maxSpecies, tOptLower, tOptUpper, compNormSd, nSpecies, repetition, speciesCombo, compMatrix,
           temperatureValues, N0, beta, delta, a_b, s, a_d, z, tmax, dt) %>% 
    rowwise() %>% 
    mutate(tOptValues = list(thermalOptFunc(tOptMin = tOptLower, tOptMax = tOptUpper, nSpecies = maxSpecies)[unlist(speciesCombo)])) %>% 
    ungroup()
  
  n <- nrow(simulationSettingsTib)
  # create empty list to store simulation results
  simulationResults <- list()
  
  # start simulation for all scenarios
  for(i in 1:nrow(simulationSettingsTib)){
    
    # create empty tibble to store simulation results
    modelResults <- as_tibble(matrix(nrow = tmax/dt + 1, ncol = simulationSettingsTib$nSpecies[i]+1))  %>% 
      mutate_if(is.logical, as.numeric) 
    
    modelResults[1,] <- t(c(0, rep(N0, times = simulationSettingsTib$nSpecies[i])))
    names(modelResults) <- c('time', paste('species', 1:simulationSettingsTib$nSpecies[i], sep = ''))
    speciesIdentity <- unlist(simulationSettingsTib$speciesCombo[i])
    
    #subset the competition matrix for simulated species
    compMatrixRun <- do.call(rbind, simulationSettingsTib$compMatrix[i])
    compMatrixRun <- compMatrixRun[speciesIdentity, speciesIdentity]
    tOptsSimulation <- round(unlist(simulationSettingsTib$tOptValues[i]),1)
    
    #subset temperatures
    temperatures <- unlist(simulationSettingsTib$temperatureValues[i])
    
    for(k in 2:nrow(modelResults)){
      
      oldBiomass <- as.numeric(modelResults[k-1, -1])  # old Biomass as vector
      newBiomass <- oldBiomass + rk4step(dt = simulationSettingsTib$dt[i], 
                                         state = round(oldBiomass,2),
                                         temp = round(temperatures[k-1], 1),
                                         b_opt = tOptsSimulation, 
                                         compMatrix = compMatrixRun,
                                         beta = beta, 
                                         delta = delta, 
                                         a_b = a_b, 
                                         s = s, 
                                         a_d = a_d, 
                                         z = z)
      
      # when biomass gets very low, species goes extinct
      newBiomass <- ifelse(newBiomass<0.01, 0, newBiomass)
      modelResults[k,] <- t(c((k-1)*simulationSettingsTib$dt[i], newBiomass))
    }
    
    modelResults <- mutate(modelResults, runID = simulationSettingsTib$runID[i])
    simulationResults[[i]] <- modelResults
    
    progress(i,n)
  }
  
  # save simulation settings and results
  simulationDat <- list(simulationSettings = simulationSettingsTib,
                        simulationResults = simulationResults)
  
  if(saveFile == T){
    saveRDS(simulationDat, paste(destPath,distType, '_R', maxSpecies,'.RData', sep = ''))
  }
  return(simulationDat)
}

# modelResults  %>% 
#   pivot_longer(values_to = 'biomass', names_to = 'species', -time) %>% 
#   ggplot(.,aes(x = time, y = biomass, colour = species)) +
#   geom_line(linewidth = 0.5) 
# 
# 


















