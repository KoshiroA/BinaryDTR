
# ====================================================================== #
# G-estimation based on collapsible measures and coherent nuisance model #
# ====================================================================== #

Gest_OP <- function(data, scale = c("add", "mult"), opt = c("max", "min"),
                    PSM = c("correct", "incorrect"), thres = 1e-6, efmA) {
  
  Ts <- max(as.numeric(gsub("A", "", grep("^A[0-9]+$", names(data), value = TRUE))))
  scale <- match.arg(scale)
  opt <- match.arg(opt)
  PSM <- match.arg(PSM)
  
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
                       "OTRAll", "Vmod","Vmod.min","Vmod.max","Vtrue", "Vsim",
                       paste0("OTRsim", Ts:1), "OTRsimAll")
  
  n <- nrow(data)
  psi_list <- list()
  cons_list <- numeric(Ts)
  H_list <- list()
  ef_list <- list()
  PS_list <- vector("list", Ts)
  
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
  
  H_list[[Ts + 1]] <- data$Y
  m_t = H_list[[Ts + 1]]
  
  for (t in Ts:1) {
    Acur <- paste0("A", t); Apre <- if (t > 1) paste0("A", t - 1) else NULL
    Xcur <- paste0("X", t); Xpre <- if (t > 1) paste0("X", t - 1) else NULL
    
    if(efmA == TRUE){
      Z_list = if(t > 1) {c(paste0("X", t), paste0("A", (t-1) ))}else if(t == 1) {paste0("X", t)}
      
    } else if(efmA == F){
      Z_list = paste0("X", t)
    }
    Z = data %>% mutate(int = 1) %>% select(int, all_of(Z_list))  %>% as.matrix()
    
    W_list = if(t > 1) {c(paste0("A", 1:(t-1)), paste0("X", 1:t ))}else if(t == 1) {paste0("X", 1:t)}
    W = data %>% mutate(int = 1) %>% select(int, W_list)  %>% as.matrix()
    H_next <- H_list[[t + 1]]
    PS_t <- PS_list[[t]]
    pars = LSEst(H_next, data[[Acur]], Z, W, scale = scale)
    psi = DREst(H_next, data[[Acur]], Z, W, pars,PS_t, scale = scale)
    
    psi_hat <- psi$x
    psi_list[[t]] <- psi_hat
    ef_t <- as.vector(Z %*% psi_hat)
    ef_list[[t]] <- ef_t
    
    d_t <- if (opt == "max") as.integer(ef_t > 0) else as.integer(ef_t < 0)
    A_t <- data[[Acur]]
    
    H_t <- if (scale == "add") {
      ifelse(d_t == 1, H_next + ef_t * (1 - A_t), H_next - ef_t * A_t)
    } else {
      ifelse(d_t == 1, H_next * exp(ef_t * (1 - A_t)), H_next * exp(-ef_t * A_t))
    }
    
    H_list[[t]] <- H_t
    
    dopt_t <- data[[paste0("d", t, ".opt")]]
    cons_list[t] <- mean(d_t == dopt_t)
    
    if(efmA == TRUE & t > 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t), paste0("A", t-1)))] <- psi_hat
    if(efmA == F | t == 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t)))] <- psi_hat
    
    if(t == 1){
      if(scale == "mult") EY = getProbRR(Z %*% psi$x, W %*% psi$alpha)
      if(scale == "add") EY = getProbRD(Z %*% psi$x, W %*% psi$alpha)
      H_list[[t]]  = ifelse(d_t == 0, EY[,1],EY[,2])
    }
  }
  
  # Final results
  res_i[1, paste0("OTR", Ts:1)] <- cons_list
  res_i[1, "OTRAll"] <- mean(Reduce(`*`, cons_list))
  res_i[1, "Vmod"] <- mean(H_list[[1]])
  res_i[1, "Vmod.min"] <- min(H_list[[1]])
  res_i[1, "Vmod.max"] <- max(H_list[[1]])
  
  return(res_i)
}

