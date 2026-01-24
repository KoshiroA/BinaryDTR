
# ====================== #
# IMPLEMENT Simulation 2 #
# ====================== #

# PARAMETERS # 
# ---------------------------------------------------------------------
load("./seed_summary.RData")
seed.summary$seed_2a.5000

n = 500; method = "M";
opt = "max"; efmA = FALSE; NLterm = FALSE
eff = c(0.3, -0.45); Ts = 4; PSM = "correct"
# ---------------------------------------------------------------------

# PREPARE for parallelization
# ---------------------------------------------------------------------
cores <- detectCores()-1
registerDoParallel(cores)
stopImplicitCluster()
# ---------------------------------------------------------------------

# ==============================================================================

for (method in c("MOP")) {
  start = proc.time()
  for (n in c(500,1000,5000)) {
    if(n == 500){set.seed(20250725)}
    if(n == 1000){seed.summary$seed_2a.500}
    if(n == 5000){seed.summary$seed_2a.1000}
    
    result <- foreach(i = 1:500, .combine = rbind, .packages = c("tidyverse","geeM","nleqslv","mgcv","brm")) %dorng% {
      
      data = datagen.cont(n = n, eff = eff, NLterm = NLterm, Ts)
      
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
        res_i = Gest_OP(data = data, opt = opt, scale = "mult", PSM = PSM, efmA = efmA)
      }
      if(method == "AOP"){
        res_i = Gest_OP(data = data, opt = opt,scale = "add", PSM = PSM, efmA = efmA)
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
    
    seed = .Random.seed
    assign(paste0("seed","_2b.",n),seed)
    rm(seed)
    atmp = paste0("res",method,"_2b.",n)
    if(PSM == "correct"){
      assign(atmp,result)
    }
    if(PSM == "incorrect"){
      assign(paste0(atmp, ".nPS"),result)
    }
    rm(atmp)
  }
  end = proc.time()
}
