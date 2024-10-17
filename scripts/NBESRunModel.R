library(sn)
source('scripts/NBESGenerateCompetitionMatrix.R')
source('scripts/NBESGenerateTOpt.R')
source('scripts/NBESGenerateTemperatures.R')
source('scripts/NBESrk4Solver.R')
source('scripts/NBESGrowthRates.R')

# solve model for 20 degrees
dt = 0.025
tmax = 200
N0 = 0.1
maxSpecies = 5

# create competition matrix
compMatrix <- generateCompetitionMatrix(competitionScenario = 'noCompetition', 
                                        nSpecies = maxSpecies,
                                        meanNormComp = 0.5,
                                        sdNormComp = 0.1,
                                        skewedDistOmegaMin = 0.2, 
                                        skewedDistOmegaMax = 0.2)

# get thermal optima for species
tOpts <- thermalOptFunc(tOptMin = 17, tOptMax = 17, nSpecies = maxSpecies)

# create temperature scenario
temperatures <- generateTemperatures(distType = 'combined', 
                                     pTempMin = 14, pTempMax = 19, 
                                     fTempMin = -1, fTempMax = 1, 
                                     #fTempMin = 14, fTempMax = 20, 
                                     fTempPeriod = 1, 
                                     tmax = tmax, dt = dt)
# temperatures <- rep(16, times = tmax/dt)
# save model settings in list:
simulationSettings <- list(runID = 1,
                           tmax = tmax,
                           dt = dt,
                           maxSpecies = maxSpecies,
                           distType = 'combined',
                           temperatures = temperatures,
                           competitionScenario = 'leftSkewedIndividualised',
                           competitionMatrix = compMatrix)


# now the simulation needs to iterate through all possible species combinations from monocultures to a 10-species community
for(i in 1:maxSpecies){
  
  nSpecies = i
  
  #all combinations of species:
  speciesComb <- t(combn(1:maxSpecies, nSpecies))
}


modelResults <- as_tibble(matrix(nrow = tmax/dt + 1, ncol = nSpecies+1))  %>% 
  mutate_if(is.logical, as.numeric) 

modelResults[1,] <- t(c(0, rep(N0, times = nSpecies)))
names(modelResults) <- c('time', paste('species', 1:nSpecies, sep = ''))
speciesIdentity <- speciesComb[1,]

#subset the competition matrix for simulated species
compMatrixSimulation <- compMatrix[speciesIdentity, speciesIdentity]
tOptsSimulation <- round(tOpts[speciesIdentity],1)

for(i in 2:nrow(modelResults)){
  
  oldBiomass <- as.numeric(modelResults[i-1, -1])  # old Biomass as vector
  newBiomass <- oldBiomass + rk4step(dt = dt, 
                                     state = round(oldBiomass,2),
                                     temp = round(temperatures[i-1], 1),
                                     b_opt = tOptsSimulation, 
                                     compMatrix = compMatrixSimulation,
                                     beta = beta, 
                                     delta = delta, 
                                     a_b = a_b, 
                                     s = s, 
                                     a_d = a_d, 
                                     z = z)
  newBiomass <- ifelse(newBiomass<0.01, 0, newBiomass)
  if(any(newBiomass/oldBiomass>5, na.rm=T)){
    break
  }
  
  modelResults[i,] <- t(c((i-1)*dt, newBiomass))
}

modelResults %>% 
  pivot_longer(values_to = 'biomass', names_to = 'species', -time) %>% 
  ggplot(.,aes(x = time, y = biomass, colour = species)) +
  geom_line(linewidth = 0.5) 


