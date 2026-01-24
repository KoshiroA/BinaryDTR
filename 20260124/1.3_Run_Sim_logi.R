# ====================== #
# IMPLEMENT Simulation 1 #
# ====================== #

# PARAMETERS # 
# ---------------------------------------------------------------------
set.seed(20250725)
load("./seed_summary.RData")

n = 500; method = "AOP"
opt = "max"; PSM = "incorrect"; NLterm = FALSE; efmA = TRUE
eff1 = -0.2; eff2 = c(0.1, 0.2, -0.1) # Sim_1a
# eff1 = 0; eff2 = c(0, 0, 0) # Sim_1n
# ---------------------------------------------------------------------

# PREPARE for the parallelization
# ---------------------------------------------------------------------
cores <- detectCores()
registerDoParallel(cores)
stopImplicitCluster()
# ---------------------------------------------------------------------

for (method in c("MN", "AN")) {
  for (n in c(500, 1000, 5000)) {
    if(n == 500){set.seed(20250725)}
    if(n == 1000){.Random.seed = seed_1a.500}
    if(n == 5000){.Random.seed = seed_1a.1000}
    
    result <- foreach(i = 1:1000, .combine = rbind, .packages = c("tidyverse","geeM","nleqslv","mgcv","brm","rlang")) %dorng% {
      
      data = datagen.logi(n = n, eff1 = eff1, eff2 = eff2, NLterm = NLterm)
      
      # implement at stage 2
      if(method == "T"){
        Vhat = datagen.cont(n = 100000)
        Vhat = Vhat$P.opt %>% mean
        return(Vhat)
      }
      if(method == "L"){
        res_i = QLearn(data, opt = "max", link = "logit", efmA = efmA)
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
        res_i = Gest(data = data, scale = "mult",opt = opt, PSM = PSM, efmA = efmA)
      }
      if(method == "A"){
        res_i = Gest(data = data, scale = "add",opt = opt, PSM = PSM, efmA = efmA)
      }
      if(method == "MN"){
        res_i = Gest(data = data, scale = "mult",opt = opt, PSM = PSM, ORM = "null", efmA = efmA)
      }
      if(method == "AN"){
        res_i = Gest(data = data, scale = "add",opt = opt, PSM = PSM, ORM = "null", efmA = efmA)
      }
      if(method == "MOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "mult", PS = PSM, efmA = efmA)
      }
      if(method == "AOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "add", PS = PSM, efmA = efmA)
      }

      res_i[1, c("Vtrue", "Vsim", paste0("OTRsim", 2:1), "OTRsimAll")] =  
        Vsim.logi(eff1 = eff1, eff2 = eff2, opt = "max", 
                  psi1 = res_i[4:5], psi2 = res_i[1:3], NLterm = NLterm) %>% apply(.,2,mean)
      
      return(res_i)
    }
    
    result = result %>% as_tibble()
    
    seed = .Random.seed
    assign(paste0("seed","_1a.",n),seed)
    rm(seed)
    atmp = paste0("res",method,"_1a.",n)
    if(PSM == "correct"){
      assign(atmp,result)
    }
    if(PSM == "incorrect"){
      assign(paste0(atmp, ".nPS"),result)
    }
    rm(atmp)
  }
}

