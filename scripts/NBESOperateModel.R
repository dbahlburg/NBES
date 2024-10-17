# housekeeping -> load packages, model function and set file paths
library(tidyverse)
library(svMisc)
library(foreach)
library(parallel)
library(progress)
source('scripts/NBESRunScenariosParallel.R')
destPath <- 'output/' # folder where simulation results will be stored
# --------------------------------------------------------------------------------------------------- #
# The following simulation settings are fixed during the model run:
# Disturbance choices:
distType <- 'press' # Disturbance type (fluctuation, press, combined)
pTempMin <- 15 # 1.1 pTempMin (minimum temperature in press disturbance)
pTempMax <- 20 # 1.2 pTempMax (maximum temperature in press disturbance)
tempControl <- 17.5 # temperature for control run
# when distType is "fluctuation" or "combined"
# fTempMin <- value (in "combined scenario", fTempMin determines the lower deviation from the press-disturbance temperature)
# fTempMax <- value (in "combined scenario", fTempMax determines the upper deviation from the press-disturbance temperature)
# fTempPeriod <- 1 (determines period of fluctuation, set to 1 day)

# community choices:
maxSpecies <- 5 # maximum species richness (5, 10)
compType <- 'normalDistribution' # (leftSkewedIndividualised, normalDistribution, noCompetition)
compNormMean <- 0 # mean value of normal distribution from which alphas will be drawn
repetitions <- 10 # number of repetitions for each combination of compNormMean and compNormSd to account for stochasticity
N0 = 0.1 # biomass of each species at t0

# parameters shaping the thermal performance curves (do not touch...(?))
beta <- 0.025 # density dependent birth rate (see Vasseur 2020)
delta <- 0.025 # density dependent death rate (see Vasseur 2020)
a_b <- 1 # temperature dependent birth rate parameter (see Vasseur 2020)
s <- 30 # temperature dependent birth rate parameter (see Vasseur 2020)
a_d <- 0.01 # temperature dependent death rate parameter (see Vasseur 2020)
z <- 0.2 # temperature dependent death rate parameter (see Vasseur 2020)

# model settings
tmax <- 150 # simulation time period
dt <- 0.05 # time step
# --------------------------------------------------------------------------------------------------- #
# the following parameters are provided as ranges, as the script will iterate through them to test different scenarios
# competition intensity:
compNormSd <- seq(0.1,0.5, length.out = 10) # sd of normal distribution from which alphas will be drawn 

# response diversity
tOptTib <- tibble(tOptScenario = paste('tOpt', 1:3, sep = ''),
                  tOptLower = c(15, 16, 17), # tibble with pairs of upper/lower temperatures over which community t_opts are distributed
                  tOptUpper = c(19, 18, 17))
# --------------------------------------------------------------------------------------------------- #
# Run model
simulatioResults <- runScenariosParallel(destPath = destPath,
             maxSpecies = maxSpecies,
             compType = compType,
             distType = distType,
             compNormMean = compNormMean,
             repetitions = repetitions,
             pTempMin = pTempMin,
             pTempMax = pTempMax,
             tempControl = tempControl,
             # fTempMin = fTempMin,
             # fTempMax = fTempMax,
             # fTempPeriod = fTempPeriod,
             N0 = N0,
             beta = beta,
             delta = delta,
             a_b = a_b,
             s = s,
             a_d = a_d,
             z = z,
             tmax = tmax,
             dt = dt,
             compNormSd = compNormSd,
             tOptTib = tOptTib,
             saveFile = T)
