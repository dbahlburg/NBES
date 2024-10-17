thermalOptFunc <- function(tOptMin,tOptMax,nSpecies){
  
  # create thermal optima equally distributed between tOptMin and tOptMax
  tOptVals <- seq(tOptMin, tOptMax, length.out = nSpecies)
  return(tOptVals)
}
