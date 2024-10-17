rk4step <- function(dt, state, temp, b_opt, compMatrix, beta, delta, a_b, s, a_d, z){
  
  k1 <- biomassChange(N = state, temp = temp, b_opt = b_opt, compMatrix = compMatrix,
                      beta = beta, delta = delta, a_b = a_b, s = s, a_d = a_d, z = z)
  k2 <- biomassChange(N = state+0.5*k1*dt, temp = temp, b_opt = b_opt, compMatrix = compMatrix,
                      beta = beta, delta = delta, a_b = a_b, s = s, a_d = a_d, z = z)
  k3 <- biomassChange(N = state+0.5*k2*dt, temp = temp, b_opt = b_opt, compMatrix = compMatrix,
                      beta = beta, delta = delta, a_b = a_b, s = s, a_d = a_d, z = z)
  k4 <- biomassChange(N = state+k3*dt, temp = temp, b_opt = b_opt, compMatrix = compMatrix,
                      beta = beta, delta = delta, a_b = a_b, s = s, a_d = a_d, z = z)
  
  return(dt * (k1 + 2*k2 + 2*k3 + k4) /6)
}
