
# ============================================================ #
# Data generation in Simulation 2: logistic-non-linear outcome #
# ============================================================ #


datagen.cont <- function(n, eff, NLterm, Ts) {
  id <- 1:n
  U <- rnorm(n, 0, 1)
  
  X <- list()
  A <- list()
  d.opt <- list()
  con <- list()
  eff.val <- list()
  reg <- list()
  
  X[[1]] <- rnorm(n, 1, 1)
  A[[1]] <- rbinom(n, 1, plogis(-0.3 + X[[1]]))
  
  for (t in 2:Ts) {
    X[[t]] <- rnorm(n, 0.5 * U + 0.6 * X[[t-1]] - 0.4 * A[[t-1]], 1)
    A[[t]] <- rbinom(n, 1, plogis(-0.3 + X[[t]] + 0.5 * A[[t-1]]))
  }
  
  if (NLterm == FALSE) {
    P.opt <- X[[1]] + 0.5 * U + 0.1
  } else {
    P.opt <- X[[1]] + 0.5 * U + 0.1 + sin(pi * X[[1]])
  }
  
  for (t in 1:Ts) {
    d.opt[[t]] <- ifelse(cbind(1, X[[t]]) %*% eff > 0, 1, 0)
    con[[t]] <- 1 - abs(d.opt[[t]] - A[[t]])
    eff.val[[t]] <- cbind(1, X[[t]]) %*% eff
    reg[[t]] <- eff.val[[t]] * (d.opt[[t]] - A[[t]])
  }
  
  pn <- P.opt - Reduce("+", reg)
  Y <- rbinom(n, 1, plogis(pn))
  
  data <- data.frame(id = id, U = U)
  for (t in 1:Ts) {
    data[[paste0("X", t)]] <- X[[t]]
    data[[paste0("A", t)]] <- A[[t]]
    data[[paste0("d", t, ".opt")]] <- d.opt[[t]]
    data[[paste0("reg", t)]] <- reg[[t]]
    data[[paste0("eff", t)]] <- eff.val[[t]]
    data[[paste0("con", t)]] <- con[[t]]
  }
  data$P.opt <- plogis(P.opt)
  data$P <- plogis(pn)
  data$Y <- Y
  
  return(data)
}

 