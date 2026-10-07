
# ===================================================== #
# Calculate the simulation-based values in Simulation 2 #
# ===================================================== #

Vsim.cont <- function(n, eff, psi_list, NLterm) {
  
  Ts <- length(psi_list)
  id <- 1:n
  U <- rnorm(n, 0, 1)
  
  dat = list()
  
  X = X_opt = A = d_opt = eff_val = cons = reg = eps = vector("list", Ts)
  
  X[[1]] = X_opt[[1]] = rnorm(n, 1, 1)
  A[[1]] <- ifelse((cbind(1, X[[1]]) %*% psi_list[[1]]) > 0, 1, 0)
  d_opt[[1]] <- ifelse((cbind(1, X[[1]]) %*% eff) > 0, 1, 0)
  eff_val[[1]] <- cbind(1, X[[1]]) %*% eff
  
  for (t in 2:Ts) {
    eps[[t]] = rnorm(n, 0, 1)
    X[[t]] <- 0.5 * U + 0.6 * X[[t-1]] - 0.4 * A[[t-1]] + eps[[t]]
    A[[t]] <- ifelse((cbind(1, X[[t]]) %*% psi_list[[t]]) > 0, 1, 0)
    
    X_opt[[t]] <- 0.5 * U + 0.6 * X_opt[[t-1]] - 0.4 * d_opt[[t-1]] + eps[[t]]
    d_opt[[t]] <- ifelse((cbind(1, X_opt[[t]]) %*% eff) > 0, 1, 0)
    eff_val[[t]] <- cbind(1, X[[t]]) %*% eff
  }
  
  if (NLterm == FALSE) {
    P.opt <- X[[1]] + 0.5 * U + 0.1
  } else {
    P.opt <- X[[1]] + 0.5 * U + 0.1 + sin(pi * X[[1]])
  }
  
  for (t in 1:Ts) {
    reg[[t]] = eff_val[[t]] * (d_opt[[t]] - A[[t]])
    cons[[t]] = 1 - abs(d_opt[[t]] - A[[t]])
  }
  
  pn <- plogis(P.opt - Reduce("+", reg))
  P.opt = plogis(P.opt)
  consAll = Reduce("*", cons)
  
  cons = bind_cols(cons)
  
  # res = tibble(P.opt, pn, A[[3]], d_opt[[3]], cons3 = cons[[3]], reg[[3]], cons2 = cons[[2]], reg[[2]], 
  #              cons1 = cons[[1]], reg[[1]], consAll)
  res = tibble(P.opt, pn, cons, consAll)
  
  return(res)
}
