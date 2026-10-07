# ====================== #
# IMPLEMENT Simulation 1 #
# ====================== #

# --------------
# Updated on 20260913

rm(list = ls())

# set your appropriate working directory via setwd()

source("./0.1_Libs.R")
source("./0.2_QLearn.R")
source("./0.3_Gestimation.R")
source("./0.4_Qk_estimation.R")
source("./0.5_Solve_Gestimation_equation.R")
source("./0.6_Gest_OP.R")
source("./1.1_Datagen_logi.R")
source("./1.2_VsimFunc_logi.R")

# PARAMETERS # 
# ---------------------------------------------------------------------
set.seed(20250725) # initial seed
load("./seed_summary.RData")

n = 500; method = "AOP"

opt = "max"; NLterm = FALSE
efmA = TRUE # whether Ak is included as a effect modifier
PSM = "incorrect"; ORM = "incorrect"

eff1 = -0.2; eff2 = c(0.1, 0.2, -0.1) # Sim_1a
# eff1 = 0; eff2 = c(0, 0, 0) # Sim_1n
# ---------------------------------------------------------------------

gc()

# PREPARE for the parallelization
# ---------------------------------------------------------------------
showConnections()
cores <- detectCores()
registerDoParallel(cores)
# ---------------------------------------------------------------------

# Before running the simulation, CHECK the seeds are appropriately set
nsim = 1000

cat("start at: ", as.character(round(Sys.time())), "\n")
for (method in c("L")) {
  
  start = proc.time()
  
  for (n in c(500, 1000, 5000)) {
    
    if(n == 500){set.seed(20250725)}
    if(n == 1000){.Random.seed = seed_1a.500}
    if(n == 5000){.Random.seed = seed_1a.1000}
    
    result <- foreach(i = 1:nsim, .combine = rbind, 
                      .packages = c("tidyverse","geeM", "geepack", "nleqslv","mgcv","brm","rlang")) %dorng% {
                         
      data = datagen.logi(n = n, eff1 = eff1, eff2 = eff2, NLterm = NLterm)
      
      # implement at stage 2
      if(method == "T"){
        Vhat = datagen.cont(n = 100000)
        Vhat = Vhat$P.opt %>% mean
        return(Vhat)
      }
      if(method == "L"){
        res_i = QLearn(data, opt = "max", link = "logit", ORM = ORM, efmA = efmA)
        # res_i = QLearn(data, opt = "max", link = "logit", efmA = efmA)
      }
      if(method == "MQ"){
        res_i = QLearn(data, opt = "max", link = "log", efmA = efmA)
      }
      if(method == "AQ"){
        res_i = QLearn(data, opt = "max", link = "identity", efmA = efmA)
      }
      if(method == "LGAM"){
        res_i = QLogi.GAM(data, opt = opt)
      }
      if(method == "M"){
        res_i = Gest(data = data, scale = "mult",opt = opt, PSM = PSM, ORM = ORM, efmA = efmA)
      }
      if(method == "A"){
        res_i = Gest(data = data, scale = "add",opt = opt, PSM = PSM, ORM = ORM, efmA = efmA)
      }
      if(method == "MN"){
        res_i = Gest(data = data, scale = "mult",opt = opt, PSM = PSM, ORM = "null", efmA = efmA)
      }
      if(method == "AN"){
        res_i = Gest(data = data, scale = "add",opt = opt, PSM = PSM, ORM = "null", efmA = efmA)
      }
      if(method == "MOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "mult", PS = PSM, ORM = ORM, efmA = efmA)
      }
      if(method == "AOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "add", PS = PSM, ORM = ORM, efmA = efmA)
      }

      res_i[1, c("Vtrue", "Vsim", paste0("OTRsim", 2:1), "OTRsimAll")] =  
        Vsim.logi(eff1 = eff1, eff2 = eff2, opt = "max", 
                  psi1 = res_i[4:5], psi2 = res_i[1:3], NLterm = NLterm) %>% apply(.,2,mean)
      
      return(res_i)
    }
    
    result = result %>% as_tibble()
    
    # if needed to update the seed values
    # seed = .Random.seed
    # assign(paste0("seed","_1a.",n),seed)
    # rm(seed)
    
    atmp = paste0("res",method,"_1a.",n)
    
    if(method != "L") {
      if(PSM == "correct" & ORM == "correct"){
        assign(atmp,result)
      }
      if(PSM == "incorrect" &  ORM == "correct"){
        assign(paste0(atmp, ".nPS"),result)
      }
      if(PSM == "correct" & ORM == "incorrect"){
        assign(paste0(atmp, ".nOR"),result)
      }
      if(PSM == "incorrect" & ORM == "incorrect"){
        assign(paste0(atmp, ".nBOTH"),result)
      }
    }else if(method == "L") {
      if(ORM == "correct"){ assign(atmp,result) }
      if(ORM == "incorrect"){ assign(paste0(atmp, ".nOR"), result) }
    }
    
    rm(atmp)
      
  }
  
    end = proc.time()
    time = end-start
    cat(method, ": ", round(time[3]/60, 2), "min; ", as.character(round(Sys.time())), "\n")
  
}

closeAllConnections()
