generateTemperatures <- function(distType, pTempMin, pTempMax, fTempMin, fTempMax, fTempPeriod, tmax, dt){
  
  # generate temperature dymamics for simulations/disturbance experiments
  if(distType == 'press'){
    tempVals <- seq(pTempMin, pTempMax, length.out = tmax/dt)
  }
  if(distType == 'fluctuation'){
    tempVals <- (fTempMax - fTempMin) * sin(2*pi/fTempPeriod * seq(dt, tmax,length.out = tmax/dt)) + (fTempMax + fTempMin)/2
  }
  if(distType == 'combined'){
    tempVals <- seq(pTempMin, pTempMax, length.out = tmax/dt) + (fTempMax - fTempMin) * sin(2*pi/fTempPeriod * seq(dt, tmax,length.out = tmax/dt)) + (fTempMax + fTempMin)/2
  }
  return(tempVals)
}