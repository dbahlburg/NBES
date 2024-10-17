generateCompetitionMatrix <- function(competitionScenario, nSpecies, meanNormComp, sdNormComp, skewedDistOmegaMin, skewedDistOmegaMax){
  if(competitionScenario == 'noCompetition'){
    compMatrix <- matrix(nrow = nSpecies, ncol = nSpecies)
    compMatrix[1:length(compMatrix)] <- 0
    diag(compMatrix) <- 1
    return(compMatrix)
  }
  if(competitionScenario == 'normalDistribution'){
    compMatrix <- matrix(nrow = nSpecies, ncol = nSpecies)
    compMatrix[1:length(compMatrix)] <- abs(rnorm(nSpecies^2, mean = meanNormComp, sd = sdNormComp))
    diag(compMatrix) <- 1
    return(compMatrix)
  }
  if(competitionScenario == 'leftSkewedIndividualised'){
    
    omegaVals <- seq(from = skewedDistOmegaMin, to = skewedDistOmegaMax, length.out = nSpecies)
    alphaValues <- abs(unlist(lapply(omegaVals, function(omega) rsn(n = nSpecies, xi = 0.05, omega = omega, alpha = 5))))
    
    compMatrix <- matrix(nrow = nSpecies, ncol = nSpecies)
    compMatrix[1:length(compMatrix)] <- alphaValues
    diag(compMatrix) <- 1
    return(compMatrix)
  }
}


