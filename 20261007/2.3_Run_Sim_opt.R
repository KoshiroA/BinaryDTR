
# ====================== #
# IMPLEMENT Simulation 2 #
# ====================== #

# 20260913: incorporate outcome regression misspecification

rm(list = ls())

# set your appropriate working directory via setwd()

source("./0.1_Libs.R")
source("./0.2_QLearn.R")
source("./0.3_Gestimation.R")
source("./0.4_Qk_estimation.R")
source("./0.5_Solve_Gestimation_equation.R")
source("./0.6_Gest_OP.R")
source("./2.1_Datagen_opt_general.R")
source("./2.2_VsimFunc_opt.R")

# PARAMETERS # 
# ---------------------------------------------------------------------
n = 500; method = "L";
opt = "max"; NLterm = FALSE
efmA = FALSE # whether Ak is included as a effect modifier
PSM = "incorrect"; ORM = "incorrect"
eff = c(0.3, -0.45) # scenario a: c(0.3, -0.45); b: c(0.1, -0.15)
Ts = 4 # N of time points
# ---------------------------------------------------------------------

# PREPARE for parallelization
# ---------------------------------------------------------------------
showConnections()
cores <- detectCores()-1
registerDoParallel(cores)
# ---------------------------------------------------------------------

# ==============================================================================

# Before running the simulation, CHECK the seeds are appropriately set

gc()
nsim = 1000; load("./seed_summary.RData")

cat("start at: ", as.character(round(Sys.time())), "\n")
for (method in c("L")) {
  
  start = proc.time()
  
  for (n in c(500,1000,5000)) {
    
    if(n == 500){.Random.seed = seed_3b.5000}
    if(n == 1000){.Random.seed = seed_3a.500}
    if(n == 5000){.Random.seed = seed_3a.1000}
    
    result <- foreach(i = 1:nsim, .combine = rbind, .errorhandling = "pass",
                      .packages = c("tidyverse","geeM", "geepack", "nleqslv", "mgcv", "brm")) %dorng% {
      
      data = datagen.cont(n = n, eff = eff, NLterm = NLterm, Ts)
      
      # implement at stage 2
      if(method == "T"){
        Vhat = datagen.cont(n = 100000)
        Vhat = Vhat$P.opt %>% mean
        return(Vhat)
      }
      if(method == "L"){
        res_i = QLearn(data, opt = "max", link = "logit", ORM = ORM,  efmA = efmA)
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
        res_i = Gest_OP(data = data, opt = opt, scale = "mult", PSM = PSM, ORM = ORM, efmA = efmA)
      }
      if(method == "AOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "add", PSM = PSM, ORM = ORM, efmA = efmA)
      }
      if(method == "MREV"){
        data$Y = 1 - data$Y
        res_i = Gest.rev(data = data,opt = opt,scale = "mult", PS = PSM)
      }
      
      psi_list = list()
      #Ts = (length(res_i) - 5)/3
      for(i in 1:Ts){
        psi_list[[i]] = res_i[str_detect(colnames(res_i), paste0("psi", i))]
      }
      
      res_i[1, c("Vtrue", "Vsim", paste0("OTRsim", Ts:1), "OTRsimAll")] =  
        Vsim.cont(n = 100000, eff = eff, psi_list = psi_list, NLterm = NLterm) %>% apply(.,2,mean)
      
      return(res_i)
    }
    
    result = result %>% as.tibble()
    
    # seed = .Random.seed
    # assign(paste0("seed","_2b.",n),seed)
    # rm(seed)
        
    atmp = paste0("res",method,"_3a.", n)
    
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
