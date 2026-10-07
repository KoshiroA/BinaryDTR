
# ====================================================================================== #
# Calculate the simulation-based values of the estimated regime in Additional Simulation #
# ====================================================================================== #

# const: a parameter controlling the event probability

Gsim_DR = function(n, res_mod, const = 1){
  
  # the estimated regime from "Gmod_DR" function
  psi2 = res_mod[str_detect(names(res_mod), "psi2")]
  psi1 = res_mod[str_detect(names(res_mod), "psi1")]
  
  # time-varying confounders follow the same distribution as "3.1_Datagen_DR.R"
  # while the treatments are assigned following the estimated regime 
  prob_L1 = plogis(0)
  L1 = rbinom(n, size = 1, prob = prob_L1)
  A1 = ifelse(cbind(1,L1) %*% psi1 >0, 1, 0)
  
  prob_L2 = plogis(-0.5 + 0.5 * L1 + 0.5 * A1 - 0.5 * (L1 * A1))
  L2 = rbinom(n, size = 1, prob = prob_L2)
  
  df_work <- data.frame(id = 1:n, L1, A1, L2)
  work_A2 = lm(id ~ A1*L1*L2,  data = df_work)
  full_mat_2 = as.vector.data.frame(model.matrix(work_A2))
  A2 = ifelse(full_mat_2 %*% psi2>0, 1, 0)
  
  prob_Y = plogis(
    const*(
      0.0 + 
        -1.5 * L1 +          
        0.5 * L2 + 
        -2.0 * (L1 * L2) +   
        
        -0.3 * A1 +
        -0.3 * A2 +
        
        0.8 * (L1 * A1) +
        0.8 * (L2 * A2) +
        0.1 * (L1 * A2) +
        -0.1 * (A1 * L2) +
        0.1 * (A1 * A2) +
        
        0.1 * (L1 * A1 * L2) +
        -0.1 * (L1 * A1 * A2) +
        0.1 * (L1 * L2 * A2) +
        -0.1 * (A1 * L2 * A2) +
        -0.1 * (L1 * A1 * L2 * A2) 
    )
    
  )
  Y = rbinom(n, size = 1, prob = prob_Y)
  
  return(mean(Y))
  
}
