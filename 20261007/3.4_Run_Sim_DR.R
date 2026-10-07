
# =============================================== #
# IMPLEMENT Additional Simulation in Section S2.4 #
# =============================================== #

rm(list = ls())

# set your appropriate working directory via setwd()

source("./0.1_Libs.R")
source("./0.2_QLearn.R")
source("./0.3_Gestimation.R")
source("./0.4_Qk_estimation.R")
source("./0.5_Solve_Gestimation_equation.R")
source("./0.6_Gest_OP.R")
source("./3.1_Datagen_DR.R")
source("./3.2_Gmod_DR.R")
source("./3.3_Gsim_DR.R")

# PREPARE for parallelization
# ---------------------------------------------------------------------
showConnections()
cores <- detectCores()-1
registerDoParallel(cores)
# ---------------------------------------------------------------------

# SPECIFY the simulation setting
# ---------------------------------------------------------------------
nsim = 1000; nobs = 5000
set.seed(20260929)
scale = "mult"
ORM = "incorrect"; PSM = "incorrect"
# ---------------------------------------------------------------------

{
start = proc.time()

result <- foreach(i = 1:nsim, .combine = rbind, 
                  .packages = c("dplyr", "tidyverse","geeM", "geepack", "nleqslv","mgcv","brm","rlang")) %dorng% {
                    
                    df_sim = datagen.DR(n = nobs, const = 1)
                    
                    res_mod = Gmod_DR(df_sim = df_sim, ORM = ORM, PSM = PSM, scale = scale)
                    
                    # bomis = res_mod
                    # bcor = res_mod
                    # bmis = res_mod
                    
                    res_sim = Gsim_DR(n=100000, res_mod, const = 1)
                    
                    if(i %% 10 == 0){
                      end = proc.time()
                      time = end-start
                      cat("i = ", i, ": ", round(time[3]/60, 2), "min; ", as.character(round(Sys.time())), "\n", 
                          file = "sim_progress.txt",append = TRUE)
                    }
                    
                    res_i = c(res_mod, Vsim = res_sim)
                    
                    return(res_i)
                  }

result = result %>% as_tibble()

end = proc.time()
time = end-start

cat("End: ", round(time[3]/60, 2), "min; ", as.character(round(Sys.time())), "\n")

}

closeAllConnections()
