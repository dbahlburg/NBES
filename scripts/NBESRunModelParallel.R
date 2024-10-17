runModelParallel <- function(nSpecies, N0, speciesCombo, compMatrix, tOptValues, temperatureValues,
                             runID, tmax, dt,
                             beta, delta, a_b, s, a_d, z){
  
  source('scripts/NBESrk4Solver.R')
  source('scripts/NBESGrowthRates.R')
  # create empty tibble to store simulation results
  modelResults <- tibble::as_tibble(matrix(nrow = tmax/dt + 1, ncol = nSpecies+1)) |>
    dplyr::mutate_if(is.logical, as.numeric) 
  
  modelResults[1,] <- t(c(0, rep(N0, times = nSpecies)))
  speciesIdentity <- unlist(speciesCombo)
  names(modelResults) <- c('time', paste('S', speciesIdentity, sep = ''))

  #subset the competition matrix for simulated species
  compMatrixRun <- do.call(rbind, compMatrix)
  compMatrixRun <- compMatrixRun[speciesIdentity, speciesIdentity]
  tOptsSimulation <- round(unlist(tOptValues),1)

  #subset temperatures
  temperatures <- unlist(temperatureValues)

  for(k in 2:nrow(modelResults)){

    oldBiomass <- as.numeric(modelResults[k-1, -1])  # old Biomass as vector
    newBiomass <- oldBiomass + rk4step(dt = dt,
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
    
    modelResults[k,] <- t(c((k-1)*dt, newBiomass))
  }

  # add runID and keep only 4 time steps per day to reduce object size
  modelResults <- dplyr::mutate(modelResults, runID = runID) |>
    dplyr::filter(time%%0.25 == 0)
  
  return(modelResults)
}














