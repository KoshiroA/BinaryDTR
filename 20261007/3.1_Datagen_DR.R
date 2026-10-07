
# ======================================================== #
# Data generation in Additional Simulation in Section S2.4 #
# ======================================================== #

# const: a parameter controlling the event probability

datagen.DR = function(n, A2 = NULL, A1 = NULL, opt = FALSE, const = 1){
  
  prob_L1 = plogis(0)
  L1 = rbinom(n, size = 1, prob = prob_L1)
  
  prob_A1 = plogis(-0.8 + 1.6 * L1)
  if(is.null(A1)) {A1 = rbinom(n, size = 1, prob = prob_A1)}
  if(opt == TRUE) { A1 = ifelse(L1 == 1, 1, 0) }
  
  prob_L2 = plogis(-0.5 + 0.5 * L1 + 0.5 * A1 - 0.5 * (L1 * A1))
  L2 = rbinom(n, size = 1, prob = prob_L2)
  
  prob_A2 = plogis(
    -0.5 + 
      0.3 * L1 + 
      0.3 * A1 + 
      0.5 * L2 + 
      0.6 * (L1 * L2) +    
      -0.8 * (L1 * A1) +    
      -0.6 * (A1 * L2) +    
      0.5 * (L1 * A1 * L2) 
  )
  
  if(is.null(A2)){A2 = rbinom(n, size = 1, prob = prob_A2)}
  if(opt == TRUE) {
    
    df_work <- data.frame(id = 1:n, L1, A1, L2) %>% 
      mutate(opt2 = case_when(
        (L1 == 0 & A1 == 0 & L2 == 1) | (L1 == 0 & A1 == 1 & L2 == 1) | 
          (L1 == 1 & A1 == 0 & L2 == 1) | (L1 == 1 & A1 == 1 & L2 == 1) ~ 1, 
        TRUE ~ 0
      ))
    A2 = df_work$opt2
  }
  
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
  
  df_sim <- data.frame(id = 1:n, L1, A1, L2, A2, Y) %>% 
    mutate(
      opt1 = case_when(L1 == 1 ~ 1, TRUE ~ 0),
      opt2 = case_when(
        (L1 == 0 & A1 == 0 & L2 == 1) | (L1 == 0 & A1 == 1 & L2 == 1) | 
          (L1 == 1 & A1 == 0 & L2 == 1) | (L1 == 1 & A1 == 1 & L2 == 1) ~ 1, 
        TRUE ~ 0
      ))
  
  return(df_sim)
}

# ---------------------------------------------------------------------------- #
# datagen.DR(n = 500000, opt = TRUE,const = 1) %>% select(Y) %>% apply(2,mean)

# datagen.DR(n = 50000) %>%
#   group_by(L1, A1, L2, A2) %>%
#   summarise(n = n()/500, mean(Y))
# ---------------------------------------------------------------------------- #

# calculate the conditional expectation of Y
# ============================================

# L1  =1; A1 = 0; L2 = 0; A2=0
# 
# plogis(
#   0.0 + 
#     -1.5 * L1 +   
#     0.5 * L2 + 
#     -2.0 * (L1 * L2) +  
#     -0.3 * A1 +
#     -0.3 * A2 +
#     0.8 * (L1 * A1) +
#     0.8 * (L2 * A2) +
#     0.1 * (L1 * A2) +
#     -0.1 * (A1 * L2) +
#     0.1 * (A1 * A2) +
#     0.1 * (L1 * A1 * L2) +
#     -0.1 * (L1 * A1 * A2) +
#     0.1 * (L1 * L2 * A2) +
#     -0.1 * (A1 * L2 * A2) +
#     -0.1 * (L1 * A1 * L2 * A2) 
# )
