
# =============================================================================== #
# Estimate dopt and calculate the algorithm-based values in Additional Simulation #
# =============================================================================== #

Gmod_DR = function(df_sim, ORM, PSM, scale){
  
  # Define correct and incorrect variable sets at stage 2
  work_Q2 = glm(Y ~ A2*A1*L1*L2, family=binomial(link="logit"), data = df_sim)
  full_mat_2 = as.vector.data.frame(model.matrix(work_Q2))
  correct_mat_2 = as.vector.data.frame(full_mat_2[,!str_detect(colnames(full_mat_2), "A2")])
  incorrect_mat_2 = correct_mat_2[, !str_detect(colnames(correct_mat_2), ":")]
  
  correct_set_2 = colnames(full_mat_2)[!str_detect(colnames(full_mat_2), "A2")][-1]
  incorrect_set_2 = correct_set_2[!str_detect(correct_set_2, ":")]
  
  correct_set.formula_2 = paste0(correct_set_2, collapse = " + ")
  incorrect_set.formula_2 = paste0(incorrect_set_2, collapse = " + ")
  
  # Propensity score model specification
  if(PSM == "correct"){exp_vals.formula = correct_set.formula_2
  }else if(PSM == "incorrect"){exp_vals.formula = incorrect_set.formula_2}
  
  PS.formula_2 = as.formula(paste0("A2 ~ ", exp_vals.formula))
  PS.mod_2 = glm(formula = PS.formula_2, data = df_sim, family = binomial(link = "logit"))
  PS_2 = predict(PS.mod_2, type = "response") 
  
  # Log OP model specification
  if(ORM == "correct"){ W2 = correct_mat_2 
  }else if(ORM == "incorrect"){ W2 = incorrect_mat_2 }
  pars_2 = LSEst(df_sim$Y, df_sim$A2, Z = correct_mat_2, W = W2, scale = scale)
  psi_2 = DREst(Y = df_sim$Y, A = df_sim$A2, Z = correct_mat_2, W = W2, 
                prepars = pars_2, PS = PS_2, scale = scale)

  # }else if(ORM == "incorrect"){
  #   work_Q2 = glm(Y ~ A2*A1*L1*L2, family=binomial(link="logit"), data = head(df_sim))
  #   work_coef = coef(work_Q2)
  #   A2_term = paste0(names(work_coef)[str_detect(names(work_coef), "A2")], collapse = " + ")
  #   mis.formula = as.formula(paste0("Y ~ ", A2_term, " + L2 + L1 + A1")) 
  #   fit_Q2 = glm(formula = mis.formula,
  #                family=binomial(link="logit"), data = df_sim)
  # }
  
  psi2 = psi_2$x
  names(psi2) = paste0("psi2_", c("int", correct_set_2))
  eff_2 = correct_mat_2 %*% psi_2$x
  opt_A2 <- ifelse(eff_2 > 0, 1, 0)
  
  V2 <- if (scale == "add") {
    ifelse(opt_A2 == 1, df_sim$Y + eff_2 * (1 - df_sim$A2), df_sim$Y - eff_2 * df_sim$A2)
  } else {
    ifelse(opt_A2 == 1, df_sim$Y * exp(eff_2 * (1 - df_sim$A2)), df_sim$Y * exp(-eff_2 * df_sim$A2))
  }
  
  # ---------------------------------------------------------
  # Stage 1 estimation
  # ---------------------------------------------------------

  # Define correct and incorrect variable sets at stage 2
  work_Q1 = glm(Y ~ A1*L1, family=binomial(link="logit"), data = df_sim)
  full_mat_1 = as.vector.data.frame(model.matrix(work_Q1))
  correct_mat_1 = as.vector.data.frame(full_mat_1[,!str_detect(colnames(full_mat_1), "A1")])
  incorrect_mat_1 = correct_mat_1[, !str_detect(colnames(correct_mat_1), "L1"), drop = FALSE] %>% as.matrix()
  
  correct_set_1 = colnames(full_mat_1)[!str_detect(colnames(full_mat_1), "A1")][-1]
  incorrect_set_1 = NULL
  
  correct_set.formula_1 = paste0(correct_set_1, collapse = " + ")
  incorrect_set.formula_1 = NULL
  
  # Propensity score model specification
  if(PSM == "correct"){exp_vals.formula = correct_set.formula_1
  PS.formula_1 = as.formula(paste0("A1 ~ ", exp_vals.formula))
  }else if(PSM == "incorrect"){
    PS.formula_1 = as.formula(paste0("A1 ~ 1" ))
    }
  
  PS.mod_1 = glm(formula = PS.formula_1, data = df_sim, family = binomial(link = "logit"))
  PS_1 = predict(PS.mod_1, type = "response") 
  
  # Log OP model specification
  if(ORM == "correct"){ W1 = correct_mat_1
  }else if(ORM == "incorrect"){ W1 = incorrect_mat_1 }
  
  pars_1 = LSEst(V2, df_sim$A1, Z = correct_mat_1, W = W1, scale = scale)
  psi_1 = DREst(Y = V2, A = df_sim$A1, Z = correct_mat_1, W = W1, 
              prepars = pars_1, PS = PS_1, scale = scale)
  
  psi1 = psi_1$x
  names(psi1) = paste0("psi1_", c("int", correct_set_1))
  eff_1 = correct_mat_1 %*% psi_1$x
  opt_A1 <- ifelse(eff_1 > 0, 1, 0)
  
  V1 <- if (scale == "add") {
    ifelse(opt_A1 == 1, V2 + eff_1 * (1 - df_sim$A1), V2 - eff_1 * df_sim$A1)
  } else {
    ifelse(opt_A1 == 1, V2 * exp(eff_1 * (1 - df_sim$A1)), V2 * exp(-eff_1 * df_sim$A1))
  }
  
  df_sim$d1 = as.vector(opt_A1)
  df_sim$d2 = as.vector(opt_A2)
  
  if(scale == "add") EY = getProbRD(correct_mat_1 %*% psi_1$x, W1 %*% psi_1$alpha)
  if(scale == "mult") EY = getProbRR(correct_mat_1 %*% psi_1$x, W1 %*% psi_1$alpha)
  df_sim$Vpred = ifelse(opt_A1 == 0, EY[,1],EY[,2])
  
  df_sim2 = df_sim %>% mutate(OTP2 = 1-(d2-opt2)^2, OTP1 = 1-(d1-opt1)^2, OTP12 = OTP2*OTP1) %>% 
    apply(2, mean)
  res = df_sim2[c("OTP2", "OTP1", "OTP12", "Y", "Vpred")]
  res = c(res, V = mean(V1), psi2, psi1)
  
  return(res)
}
