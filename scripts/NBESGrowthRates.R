# temperature-dependent rate of biomass change
biomassChange <- function(N, b_opt, temp, compMatrix, beta, delta, a_b, s, a_d, z){

  # b0T: temperature dependent birth rate (Vasseur 2020)
  b0T = a_b * exp(-(temp - b_opt)^2/s)
  
  # d0T: temperature-dependent death rate (Vasseur 2020)
  d0T <- a_d * exp(z*temp)
  
  # derive temperature-dependent growth rates r(T) and capacity k(T) to simplify growth equation
  rT <- b0T - d0T
  kT <- (b0T - d0T)/(beta+delta)
  
  if(length(compMatrix)==1){
    #calculate growth rate
    return(rT*N*(1 - N/kT))
  } else {
    return(rT*N*(1 - colSums(compMatrix * N/kT)))
  }
  
}
