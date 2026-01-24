
# ======================================================== #
# Data generation in Simulation 1: logistic-linear outcome #
# ======================================================== #

# rm(list = ls())

datagen.logi = function(n, eff1, eff2, opt = "max", NLterm){
  id = 1:n
  U = rnorm(n, 0, 1)
  X1 = rnorm(n, 1, 1)
  A1 = rbinom(n,1,expit(-0.3 + X1))
  X2 = rnorm(n, 0.5*U + 0.6*X1 - 0.4*A1, 1)
  A2 = rbinom(n,1,expit(-0.3 + X2 + 0.5*A1))
  
  Z2 = cbind(1,X2,A1)
  Z1 = if(length(eff1) == 1) cbind(A1) else cbind(1,X1)
  theta = c(0.25, 0.5, 0.5, eff1, 0, eff2)
  core = cbind(1, U, X1, Z1*A1, X2, A2, X2*A2, A1*A2)
  
  if(NLterm == T){
    theta = c(theta, -0.1, 0.1)
    phi1 = X1^2; phi2 = X2^2
    core = cbind(core, phi1, phi2)
  }

  P = expit(core %*% theta)
  Y = rbinom(n,1,P)
  
  if(opt == "max"){
    d2.opt = ifelse(Z2 %*% eff2 > 0,1,0)
    d1.opt = ifelse(Z1 %*% eff1 > 0,1,0)
  }
  
  if(opt == "min"){
    d2.opt = ifelse(Z2 %*% eff2 < 0,1,0)
    d1.opt = ifelse(Z1 %*% eff1 > 0,1,0)
  }
  
  con1 = 1 - abs(d1.opt- A1)
  con2 = 1 - abs(d2.opt- A2)
  
  data = tibble(U, X1, A1, X2, A2, Y, d2.opt, d1.opt, P, id, con1, con2)
  
  return(data)
  
}

# CHECK the dataset
# ------------------------------------------------------------------------------
# data = datagen.logi(n = 100000, eff1, eff2, NLterm = F, opt = "max")
# mean(data$Y)
# 
# # 1. There is no extreme trend of the optimal treatment 
# # This is essentially AT or NT regime
# data %>% group_by(A1,d1.opt) %>% summarise(n())
# data %>% group_by(A2,d2.opt) %>% summarise(n())
# 
# # 2. There are inconsistencies between observed and optimal treatments
# data %>% group_by(con1,con2) %>% summarise(n())
# 
# # 3. Once treat, continue to treat
# # There are few cases of ceasing to (un)treat
# data %>% group_by(d1.opt,d2.opt) %>% summarise(n())
# 
# data %>% group_by(A1,A2) %>% summarise(n())
# 
# # 5. moderate correlation of the covariates or trt effects across time points
# pairs.panels(data[,c(2,4)],
#              method = "pearson",  # 相関係数の種類
#              hist.col = "lightblue",
#              density = TRUE,
#              ellipses = TRUE)
# ------------------------------------------------------------------------------

# TRUE value of "Always Treat" or "Never Treat" 
# ------------------------------------------------------------------------------
# gcomp.logi = function(n, eff1, eff2, opt, NLterm){
#   res_list = list()
#   
#   for (regime in c("AT","NT","OPT")) {
#     id = 1:n
#     U = rnorm(n, 0, 1)
#     X1 = rnorm(n, 1, 1)
#     #X1 = rbinom(n, 1, 0.5)
#     Z1 = if (length(eff1) == 1) cbind(rep(1,n)) else cbind(1,X1)
#     A1 = if (regime == "AT") 1 else if (regime == "NT") 0 else ifelse(Z1 %*% eff1 > 0,1,0)
#     X2 = rnorm(n, 0.5*U + 0.6*X1 - 0.4*A1, 1)
#     Z2 = cbind(1,X2,A1)
#     A2 = if (regime == "AT") 1 else if (regime == "NT") 0 else ifelse(Z2 %*% eff2 > 0,1,0)
#     
#     theta = c(0.25, 0.5, 0.5, eff1, 0, eff2)
#     core = cbind(1, U, X1, Z1*as.numeric(A1), X2, A2, X2*A2, A1*A2)
#     
#     if(NLterm == T){
#       theta = c(theta, -0.1, 0.1)
#       phi1 = X1^2; phi2 = X2^2
#       core = cbind(core, phi1, phi2)
#     }
#     
#     res_list[[regime]] = expit(core %*% theta) %>% mean()
#   }
#   return(res_list)
# }
# ------------------------------------------------------------------------------