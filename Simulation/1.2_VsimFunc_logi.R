
# ===================================================== #
# Calculate the simulation-based values in Simulation 1 #
# ===================================================== #

# TRUE coef1 and coef2 provide TRUE value function

# INPUT
# --------------------------
# eff: TRUE parameters
# psi: ESTIMATED parameters
# --------------------------

Vsim.logi = function(eff1, eff2, opt, psi1, psi2, NLterm){
  
  psi1 = as.numeric(psi1); psi2 = as.numeric(psi2)
  
  n  = 100000
  id = 1:n
  U = rnorm(n, 0, 1)
  X1 = rnorm(n, 1, 1)
  A1 = ifelse(cbind(1,X1) %*% psi1 > 0,1,0)
  lnpe1 = cbind(1,X1) %*% psi1
  d1.opt = ifelse(cbind(1,X1) %*% c(eff1, 0) > 0,1,0)
  
  # GENERATE the shared error
  eps2 = rnorm(n, 0, 1)
  
  # ESTIMATED REGIME
  X2 = 0.5*U + 0.6*X1 - 0.4*A1 + eps2
  A2 = ifelse(cbind(1,X2,A1) %*% psi2 > 0,1,0)
  
  # TRUE REGIME
  X2.opt = 0.5*U + 0.6*X1 - 0.4*d1.opt + eps2
  d2.opt = ifelse(cbind(1,X2.opt,d1.opt) %*% eff2 > 0,1,0)
  
  # FOLLOW the latter course to estimate E[Y(d^opt)]
  Z2 = cbind(1,X2,A1)
  Z1 = if(length(eff1) == 1) cbind(A1) else cbind(1,X1)
  theta = c(0.25, 0.5, 0.5, eff1, 0, eff2)
  core = cbind(1, U, X1, A1, X2, A2, X2*A2, A1*A2)
  core.opt = cbind(1, U, X1, d1.opt, X2.opt, d2.opt, X2.opt*d2.opt, d1.opt*d2.opt)
  
  lnpe2 = Z2 %*% psi2
  lnpt2 = cbind(1,X2.opt,d1.opt) %*% eff2
  
  if(NLterm == T){
    theta = c(theta, -0.1, 0.1)
    phi1 = X1^2; phi2 = X2^2
    core = cbind(core, phi1, phi2)
    core.opt = cbind(core.opt, phi1, phi2)
  }
  P = expit(core %*% theta)
  P.opt = expit(core.opt %*% theta)
  
  # CALCULATE the consistency of the two regime
  con1 = 1 - abs(d1.opt- A1)
  con2 = 1 - abs(d2.opt- A2)
  
  res = tibble(P.opt, P, con2, con1) %>%
    mutate(con12 = ifelse(con1 == 1 & con2 == 1, 1, 0))
  
  return(res)
}