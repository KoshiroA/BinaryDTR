
# ========================================================================= #
# G-estimation based on collapsible measures with GLM treatment-free models #
# ========================================================================= #

# 20260913
# Incorporate ORM = c("null", "correct", "incorrect")

Gest <- function(data, scale = c("add", "mult"), opt = c("max", "min"),
                 PSM = c("correct", "incorrect"), ORM = c("null", "correct", "incorrect"), 
                 thres = 1e-6, efmA) {
  
  Ts <- max(as.numeric(gsub("A", "", grep("^A[0-9]+$", names(data), value = TRUE))))
  scale <- match.arg(scale)
  opt <- match.arg(opt)
  PSM <- match.arg(PSM)
  
  # To store the result
  if (efmA == TRUE) {
    ncol = Ts*5 + 7 
    coef_names <- unlist(lapply(Ts:1, function(t) paste0("psi", t, ".", c("0", paste0("X", t), paste0("A", t-1)))))
  } else if (efmA == F) {
    ncol = Ts*4 + 7
    coef_names <- unlist(lapply(Ts:1, function(t) paste0("psi", t, ".", c("0", paste0("X", t)))))
  }
  res_i <- matrix(NA, nrow = 1, ncol = ncol)
  colnames(res_i) <- c(coef_names,
                       paste0("OTR", Ts:1),
                       "OTRAll", "Vmod","Vmod.min","Vmod.max", "Vtrue", "Vsim",
                       paste0("OTRsim", Ts:1), "OTRsimAll")
  
  n <- nrow(data)
  psi_list <- list()
  cons_list <- numeric(Ts)
  H_list <- list()
  V_list <- list()
  ef_list <- list()
  PS_list <- vector("list", Ts)
  
  # PS model fit and estimation across stages
  # --------------------------------------------------------------------
  for (t in Ts:1) {
    Acur <- paste0("A", t)
    Xcur <- paste0("X", t)
    Apre <- if (t > 1) paste0("A", t - 1) else NULL
    Xpre <- if (t > 1) paste0("X", t - 1) else NULL
    
    if(efmA == TRUE){
      if (t > 1) Z <- as.matrix(cbind(1, data[[Xcur]], data[[Apre]]))
      if (t == 1) Z <- as.matrix(cbind(1, data[[Xcur]]))
    } else if(efmA == F){
      Z <- as.matrix(cbind(1, data[[Xcur]]))
    }
    
    if (PSM == "correct") {
      formula <- as.formula(paste0(Acur, "~", Xcur,
                                   if (t > 1) {paste0(" + ", paste0("A", 1:(t-1), collapse = " + "),
                                                      " + ", paste0("X", 1:(t-1), collapse = " + "))} else ""))
    } else {
      formula <- as.formula(paste0(Acur, "~", 
                                   if (t > 1) {paste0("A", 1:(t-1), collapse = " + ")} else "1"))
      
    }
    
    PS_model <- glm(formula, data = data, family = binomial())
    PS_list[[t]] <- predict(PS_model, type = "response")
  }
  # --------------------------------------------------------------------
  
  # Iterative outcome model fit and estimation across stages
  # --------------------------------------------------------------------
  
  # Initialize the outcome
  H_list[[Ts + 1]] <- data$Y
  V_list[[Ts + 1]] <- data$Y
  V_t = V_list[[Ts + 1]]
  
  for (t in Ts:1) {
    Acur <- paste0("A", t); Apre <- if (t > 1) paste0("A", t - 1) else NULL
    Xcur <- paste0("X", t); Xpre <- if (t > 1) paste0("X", t - 1) else NULL
    
    # Z: variables as effect modifiers
    if(efmA == TRUE){
      Z_list = if(t > 1) {c(paste0("X", t), paste0("A", (t-1) ))}else if(t == 1) {paste0("X", t)}
      
    } else if(efmA == F){
      Z_list = paste0("X", t)
    }
    Z = data %>% mutate(int = 1) %>% select(int, all_of(Z_list))  %>% as.matrix()
    
    # W: variables included in nuisance outcome models
    
    if(ORM == "correct"){
      W_list = if(t > 1) {c(paste0("A", 1:(t-1)), paste0("X", 1:t ))}else if(t == 1) {paste0("X", 1:t)}
      }
    if(ORM == "incorrect"){
      W_list = if(t > 1) {c(paste0("A", 1:(t-1)))}else if(t == 1) {NULL}
      }
    
    W = data %>% mutate(int = 1) %>% select(int, W_list)  %>% as.matrix()
    
    # combined the two colmns
    ZW = cbind(W, data[[Acur]]*Z)
    
    H_next <- H_list[[t + 1]]
    V_next <- V_list[[t + 1]]
    PS_t <- PS_list[[t]]
    
    # Nuisance model estimation via least square
    rss_fn <- function(beta, X, y) {
      mu <- expit(X %*% beta)
      sum((y - mu)^2)
    }
    
    # Minimize the residual sum of square
    mt_opt <- optim(par = rep(0, ncol(ZW)), fn = rss_fn, X = ZW, y = H_next)
    
    ZW = cbind(W, 0*Z)
    m_t <- expit(ZW %*% mt_opt$par)
    
    if(ORM == "null"){m_t = 0}
    
    # G-estimating equation
    gest_fun <- if (scale == "add") {
      function(psi) {
        t(Z) %*% as.vector((data[[Acur]] - PS_t) * (H_next - data[[Acur]] * (Z %*% psi) - m_t))
      }
    } else {
      function(psi) {
        t(Z) %*% as.vector((data[[Acur]] - PS_t) * (H_next * exp(-data[[Acur]] * (Z %*% psi)) - m_t))
      }
    }
    
    # solve
    psi_start <- rep(0, ncol(Z))
    opt_result <- nleqslv(psi_start, gest_fun, method = "Newton", control = list(allowSingular = TRUE))
    psi_hat <- opt_result$x
    psi_list[[t]] <- psi_hat
    ef_t <- as.vector(Z %*% psi_hat)
    ef_list[[t]] <- ef_t
    
    d_t <- if (opt == "max") as.integer(ef_t > 0) else as.integer(ef_t < 0)
    A_t <- data[[Acur]]
    
    if (scale == "add") {
      H_t <- ifelse(d_t == 1, H_next + ef_t * (1 - A_t), H_next - ef_t * A_t)
      V_t <- m_t + d_t * ef_t
    } else {
      H_t <- ifelse(d_t == 1, H_next * exp(ef_t * (1 - A_t)), H_next * exp(-ef_t * A_t))
      V_t <- m_t * exp(d_t * ef_t)
    }
    
    V_list[[t]] <- V_t
    H_list[[t]] <- H_t
    
    dopt_t <- data[[paste0("d", t, ".opt")]]
    cons_list[t] <- mean(d_t == dopt_t)
    
    if(efmA == TRUE & t > 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t), paste0("A", t-1)))] <- psi_hat
    if(efmA == F | t == 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t)))] <- psi_hat
    
  }
  
  # --------------------------------------------------------------------
  
  # Final result
  res_i[1, paste0("OTR", Ts:1)] <- cons_list
  res_i[1, "OTRAll"] <- mean(Reduce(`*`, cons_list))
  res_i[1, "Vmod"] <- mean(H_list[[1]])
  res_i[1, "Vmod.min"] <- min(H_list[[1]])
  res_i[1, "Vmod.max"] <- max(H_list[[1]])
  
  return(res_i)
}

