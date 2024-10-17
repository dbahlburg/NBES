runScenariosParallel <- function(destPath, maxSpecies, distType, compType, pTempMin = NA, pTempMax = NA, tempControl, fTempMin = NA, fTempMax = NA, fTempPeriod = NA, compNormMean, repetitions, N0, beta, delta, a_b, s, a_d, z, tmax, dt, compNormSd, tOptTib, saveFile = T){
  
  # load auxiliary functions to run the model
  source('scripts/NBESGenerateCompetitionMatrix.R')
  source('scripts/NBESGenerateTOpt.R')
  source('scripts/NBESGenerateTemperatures.R')
  
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
                                 comSdId = paste('sd',k, sep = ''),
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
    crossing(speciesCombTib)  %>% 
    rowwise() %>% 
    mutate(speciesCombiMerged = paste('S', unlist(speciesCombo), sep='', collapse = '')) %>% 
    ungroup() %>% 
    mutate(#runID = row_number(),
           runID = paste(comSdId, tOptScenario, paste('R', repetition,sep=''), speciesCombiMerged, sep = '_'),
           communityID = paste(comSdId, tOptScenario, paste('R', repetition,sep=''), sep = '_')) %>% 
    crossing(tempTib) %>% 
    mutate(runID = paste(runID, distTreatment, sep='_'),
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
    select(runID, communityID, distType, distTreatment, compType, compNormMean, maxSpecies, tOptLower, tOptUpper, compNormSd, nSpecies, repetition, speciesCombo, compMatrix,
           temperatureValues, N0, beta, delta, a_b, s, a_d, z, tmax, dt) %>% 
    rowwise() %>% 
    mutate(tOptValues = list(thermalOptFunc(tOptMin = tOptLower, tOptMax = tOptUpper, nSpecies = maxSpecies)[unlist(speciesCombo)])) %>% 
    ungroup()
  
  iterations <- nrow(simulationSettingsTib)
  ## ------------------------------------------------------------------------------------- #
  # set up "infrastructure" to parallelize processing
  # detect number of cores (subtract 1 for safety)
  nCores <- parallel::detectCores()
  
  #create the cluster
  myCluster <- parallel::makeCluster(
    nCores, 
    type = "PSOCK"
  )
  
  #register cluster to be used by %dopar%
  doSNOW::registerDoSNOW(cl = myCluster)
  # ------------------------------------------------------------------------------------- #
  # create progress bar 
  pb <- progress_bar$new(
    format = "[:bar] :current/:total | :elapsed | eta: :eta | at :tick_rate/sec",
    total = iterations,
    show_after = 0,
    width = 60)
  # allowing progress bar to be used in foreach -----------------------------
  progress <- function(n){
    pb$tick()
  } 
  
  opts <- list(progress = progress)
  # ------------------------------------------------------------------------------------- #
  # start the parallelized loop
  simulationResults <- foreach(i = 1:iterations, .options.snow = opts) %dopar% {
  
    # load function to run model parallel
    source('scripts/NBESRunModelParallel.R')
    runModelParallel(N0 = simulationSettingsTib$N0[i],
                     tmax = simulationSettingsTib$tmax[i],
                     dt = simulationSettingsTib$dt[i],
                     nSpecies = simulationSettingsTib$nSpecies[i],
                     speciesCombo = simulationSettingsTib$speciesCombo[i],
                     compMatrix = simulationSettingsTib$compMatrix[i],
                     tOptValues = simulationSettingsTib$tOptValues[i],
                     temperatureValues = simulationSettingsTib$temperatureValues[i],
                     runID = simulationSettingsTib$runID[i],
                     beta = simulationSettingsTib$beta[i], 
                     delta = simulationSettingsTib$delta[i], 
                     a_b = simulationSettingsTib$a_b[i], 
                     s = simulationSettingsTib$s[i], 
                     a_d = simulationSettingsTib$a_d[i], 
                     z = simulationSettingsTib$z[i])
  }
  
  # close the cluster
  parallel::stopCluster(cl = myCluster)
  
  # save simulation settings and results
  simulationDat <- list(simulationSettings = simulationSettingsTib,
                        simulationResults = simulationResults)
  
  if(saveFile == T){
    saveRDS(simulationDat, paste(destPath,distType, '_R', maxSpecies,'.RData', sep = ''))
  }
  return(simulationDat)
}
