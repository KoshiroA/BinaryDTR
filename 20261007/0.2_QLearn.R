
# ======================================= #
# Q-learning with standard link functions #
# ======================================= #

# 20260913
# Incorporate ORM = c("correct", "incorrect")

QLearn <- function(data, opt = c("max", "min"), link = c("logit", "log", "identity"), 
                   efmA, ORM = c("correct", "incorrect")) {
  
  # OBTAIN the number of time points
  Ts <- max(as.numeric(gsub("A", "", grep("^A[0-9]+$", names(data), value = TRUE))))
  
  opt <- match.arg(opt)
  link <- match.arg(link)
  expit <- function(x) 1 / (1 + exp(-x))
  logit <- function(p) log(p / (1 - p))
  
  link = "logit"
  inv_link <- switch(link,
                     logit = function(x) 1 / (1 + exp(-x)),
                     log = function(x) exp(x),
                     identity = function(x) x
  )
  family_used = 
    if (link == "identity") {gaussian(link = link)} else if (link == "log") {poisson(link = link)} else {binomial(link = link)}
  
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
  
  H <- data$Y
  decision_list <- list()
  coef_list <- list()
  
  for (t in Ts:1) {
    Xcur <- paste0("X", t)
    Acur <- paste0("A", t); Apre <- ifelse(t>1, paste0("A", t-1), NA)
    
    # Z: variables as effect modifiers
    if(efmA == TRUE){
      if (t > 1) {
        Z <- as.matrix(cbind(1, data[[Xcur]], data[[Apre]]))
        Z_nam = c(Xcur, Apre)
      }
      if (t == 1) {
        Z <- as.matrix(cbind(1, data[[Xcur]]))
        Z_nam = c(Xcur)
      }

    } else if(efmA == F){
      Z <- as.matrix(cbind(1, data[[Xcur]]))
      Z_nam = c(Xcur)
    }
    efm = paste0(Acur, ":", Z_nam, collapse = " + ")
    
    if(ORM == "correct"){
      main = paste0(
        paste0("A", t:1, collapse = " + "), " + ", paste0("X", t:1, collapse = " + ")
      )
    }
    if(ORM == "incorrect"){
      main = paste0("A", t:1, collapse = " + ")
    }
    
    formula <- as.formula(paste0("H ~ ", efm," + ", main))
    
    # geeglm or geem
    if (t == Ts) {
      # formula <- as.formula(paste0("H ~ ", efm," + ", main))
  
                                   # #if (Ts > 1) paste0(if (efmA == TRUE) paste0("+", Acur, "*", Apre),
                                   #                   " + ", paste0("A", 1:(Ts-1), collapse = " + "),
                                   #                   " + ", paste0("X", 1:(Ts-1), collapse = " + ")) else ""))
      mod <- glm(formula, family = family_used, data = data)
      coef_all <- coef(mod)
      # coef_int = if (efmA == F) coef_all[c(2, length(coef_all))] else coef_all[c(2, length(coef_all)-1,length(coef_all))]
      
      coef_int = coef_all[c(Acur, paste0(Acur, ":", Z_nam))]
      
    } else {
      # formula <- as.formula(paste0("H ~ ", Acur, "*", Xcur,
      #                              if (t > 1) paste0(if (efmA == TRUE) paste0("+", Acur, "*", Apre),
      #                                                " + ", paste0("A", 1:(t-1), collapse = " + "),
      #                                                " + ", paste0("X", 1:(t-1), collapse = " + ")) else ""))
      
      mod <- geem(formula, family = family_used, data = data, id = data[["id"]])
      
      
      coef_all <- coef(mod)
      
      # if (t == 1) coef_int = if (efmA == F) coef_all[c(2, length(coef_all))] else coef_all[c(2,length(coef_all))]
      # if (t > 1) coef_int = if (efmA == F) coef_all[c(2, length(coef_all))] else coef_all[c(2, length(coef_all)-1,length(coef_all))]
      
      coef_int = coef_all[c(Acur, paste0(Acur, ":", Z_nam))]
    }
    
    ef <- Z %*% coef_int
    d <- if (opt == "max") ifelse(ef > 0, 1, 0) else ifelse(ef < 0, 1, 0)
    
    if (t == Ts) {
      
      pred.data = data %>% select(paste0(c("A","X"), Ts), paste0("A", Ts:1), paste0("X", Ts:1)) %>% 
        mutate(!!names(.)[1] := d)
      
      Xmat <- model.matrix(formula, data = pred.data)
      H_new <- as.vector(inv_link(Xmat %*% coef_all))
      
      H_new <- if (link == "identity") {
        ifelse(d == 1, H + ef * (1 - data[[Acur]]), H - ef * data[[Acur]])
      } else if (link == "log"){
        ifelse(d == 1, H * exp(ef * (1 - data[[Acur]])), H * exp(-ef * data[[Acur]]))
      } else if (link == "logit"){as.vector(inv_link(Xmat %*% coef_all))}
      
      cons <- 1 - (data[[paste0("d", t, ".opt")]] - d)^2
    } else {
      ef_term <- rowSums(sapply((t+1):Ts, function(k) {
        if (!is.null(decision_list[[as.character(k)]]))
          decision_list[[as.character(k)]]$ef * decision_list[[as.character(k)]]$d
        else 0
      }))
      
      pred.data = data %>% select(paste0(c("A","X"), t),paste0("A", t:1),paste0("X", t:1)) %>% 
        mutate(!!names(.)[1] := d)
      Xmat <- model.matrix(formula, data = pred.data)
      
      # H_new <- as.vector(inv_link(as.matrix(cbind(1, pred.data)) %*% as.vector(mod$beta)))
      
      H_new <- if (link == "identity") {
        ifelse(d == 1, H_new + ef * (1 - data[[Acur]]), H_new - ef * data[[Acur]])
      } else if (link == "log"){
        ifelse(d == 1, H_new * exp(ef * (1 - data[[Acur]])), H_new * exp(-ef * data[[Acur]]))
      } else if (link == "logit"){as.vector(inv_link(Xmat %*% coef_all))}
      
      # H_new <- 
      cons <- 1 - (data[[paste0("d", t, ".opt")]] - d)^2
    }
    
    data[[paste0("H", t)]] <- H_new
    data[[paste0("d", t)]] <- d
    data[[paste0("ef", t)]] <- ef
    data[[paste0("cons", t)]] <- cons
    H <- H_new
    
    decision_list[[as.character(t)]] <- list(ef = ef, d = d)
    coef_list[[as.character(t)]] <- coef_int
    
    # save the result
    if(efmA == TRUE & t > 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t), paste0("A", t-1)))] <- coef_int
    if(efmA == F | t == 1) res_i[1, paste0("psi", t, ".", c("0", paste0("X", t)))] <- coef_int
    res_i[1, paste0("OTR", t)] <- mean(cons)
  }
  
  data$consAll <- rowSums(sapply(1:Ts, function(t) data[[paste0("cons", t)]])) == T
  res_i[1, "OTRAll"] <- mean(data$consAll)
  res_i[1, "Vmod"] <- mean(data[[paste0("H1")]])
  res_i[1, "Vmod.min"] <- min(data[[paste0("H1")]])
  res_i[1, "Vmod.max"] <- max(data[[paste0("H1")]])
  
  return(res_i)
}

