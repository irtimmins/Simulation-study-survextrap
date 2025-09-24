#############################################################
# Load packages.
#############################################################

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(survextrap)
library(flexsurv)
library(rslurm)
library(rsimsum)
library(stringr)
library(readr)
library(rstan)
library(tibble)
library(purrr)
library(forcats)
library(data.table)
library(emg)
library(rstpm2)

#############################################################
# specify jobname and path where results are to be stored.
#############################################################

# set to project root.
# setwd("/home/klvq491/Documents/simulation_survextrap/simsurvextrap/")

jobname <- "trt5_rw" 
user <- Sys.info()["user"]
scratch_directory <- paste0("/wscratch/", user, "/")

store_res <- paste0(scratch_directory, "simsurvextrap_slurm_", jobname, "/") 



#############################################################
# import trial data and call scripts to create functions.
#############################################################

source("R/simulate_dgm_treatment_effect.R")
source("R/simulate_dgm.R")
source("R/functions_treatment_effect.R")
source("R/fit_model.R")
source("R/estimands.R")
source("R/visualise.R")
source("R/performance.R")

setwd(store_res)


#######################################################
# Specify survextrap scenarios
########################################################

# write scenarios as tibble

nsim <- 1000 # simulation replicates


scenarios_flex <- expand_grid(dgm = 1:3,
                              k_knots = c(0,1,2,4,8)) %>%
  mutate(flex_scenario_id = row_number()) %>%
  relocate(flex_scenario_id)

for(i in 1:nrow(scenarios_flex)){
  if(!dir.exists(paste0("scen_flex", i))){
    dir.create(paste0("scen_flex", i))
  }  
}

#######################################################
# create tibble/data.frame of parameters to call from
#######################################################

set.seed(132497059)

pars_flex <- expand_grid(flex_scenario_id = scenarios_flex$flex_scenario_id,
                         isim = seq(nsim)) %>%
  left_join(scenarios_flex, by = join_by(flex_scenario_id)) %>%
  mutate(data_file = paste0("../dgm", dgm, "/sim_data", isim ,".rds")) %>%
  mutate(save_file = paste0("../scen_flex", flex_scenario_id, "/flex_est", isim,".rds")) %>%
  mutate(seed =  round(runif(n(), min = 10*7, max = 10^9))) %>%
  select(-c(flex_scenario_id, dgm, isim)) %>%
  slice(sample(1:n())) 
saveRDS(pars_flex, file = paste0("pars_flex_",jobname, ".rds"))

head(pars_flex)
#length(unique(pars_flex$save_file))
#length(unique(pars_flex$data_file))

# test flex_fit

#temp <- readRDS("scen_flex1/flex_est1.rds",)
if(!dir.exists("test")){
  dir.create("test")
}  

flex_fit(k_knots = 4, 
         data_file ="dgm2/sim_data10.rds",
         save_file = "test/flex_test.rds", 
         seed = 100)

test <- readRDS("test/flex_test.rds")
test
#temp
#temp
test <- readRDS("_rslurm_flex/results_14.RDS")
test

flex_fit <- function(k_knots, data_file, save_file, seed){
  
  set.seed(seed) 
  
  data <- readRDS(data_file)
  
  data <- data %>%
    mutate(trt = as.factor(trt))
  
  if(k_knots == 0) anc_list <- list(gamma1=~trt)
  
  if(k_knots == 1) anc_list <- list(gamma1=~trt, gamma2=~trt)
  
  if(k_knots == 2) anc_list <- list(gamma1=~trt, gamma2=~trt,
                                    gamma3=~trt)
  
  if(k_knots == 4) anc_list <- list(gamma1=~trt, gamma2=~trt,
                                    gamma3=~trt, gamma4=~trt, gamma5=~trt)
  
  if(k_knots == 8) anc_list <- list(gamma1=~trt, gamma2=~trt,
                                    gamma3=~trt, gamma4=~trt, gamma5=~trt,
                                    gamma6=~trt, gamma7=~trt, gamma8=~trt, 
                                    gamma9=~trt)
  
  flexmod0 <- try(flexsurvspline(Surv(time, event)~trt, 
                                 data = data , 
                                 k = k_knots,
                                 anc=anc_list,
                                 spline = "splines2ns"))
  
  mod0_run_status <- inherits(flexmod0, "try-error")
  
  
  # model run succesfully, 
  # continue with standsurv:
  
  if(!mod0_run_status){
    
    covar_status <- !is.na(flexmod0$cov[1])
    
    if(covar_status){
      
      est_irmst <- standsurv(flexmod0,at = list(list(trt = as.factor(0)),list(trt = as.factor(1))), 
                             type = "rmst", t = 5, contrast = "difference", ci = TRUE, se = TRUE, boot = T, B = 1000) %>%
        mutate(estimand = "irmst", trt = NA,
               value = contrast2_1, 
               value_ci_low = contrast2_1_lci, 
               value_ci_high = contrast2_1_uci,
               value_se = contrast2_1_se,
               t= time) %>%
        select(estimand, trt, t, value, value_ci_low, value_ci_high, value_se)
      
      
      est <- est_irmst
      
      saveRDS(est, save_file)
      
    }
  }
}

head(pars_flex)

sjob_flex <- slurm_apply(flex_fit , 
                         pars_flex, 
                         jobname ="flex_trt",
                         nodes = 40, 
                         cpus_per_node = 4, 
                         submit = TRUE,
                         pkgs = c("dplyr", "tidyr","readr",
                                  "flexsurv", "survextrap", "rstan","rstpm2" ),
                         slurm_options = list(time='04:00:00',
                                              partition='core',
                                              "mem-per-cpu"= '4G'))

check_status <- paste0("sacct -S ", as.character(Sys.Date()-4),
                       " -u ", user ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")
system(check_status)

temp <- readRDS("_rslurm_flex_trt/results_0.RDS")
temp


for(i in scenarios_flex$flex_scenario_id){
  est <- NULL
  
  for(j in 1:nsim){
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
    if(file.exists(paste0("scen_flex", i, "/flex_est", j,".rds"))){
      temp <- readRDS(paste0("scen_flex", i, "/flex_est", j,".rds"))  %>%
        filter(estimand == "irmst") %>%
        mutate(isim = j, mod_type = "sim") 
      
      
      est <- rbind(est, temp)
    }
    
    
  }
  #nrow(est)
  
  saveRDS(est, paste0("scen_flex", i, "_est_rmst.rds"))
  
  true <- readRDS(paste0("dgm",scenarios_flex$dgm[i],"/true_est.RDS"))
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true") %>%
    filter(t == 5) %>%
    filter(estimand == "irmst")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("scen_flex", i, "_res_rmst.rds"))
  
}

#head(res)






####################################################
# benchmarking for rstpm2
####################################################


# write scenarios as tibble

nsim <- 1000 # simulation replicates

scenarios_rstpm <- expand_grid(dgm = 1:3,
                               df_base = c(1,2,3,5,9),
                               df_tvc = c(1,2,3,5,9)) %>%
  mutate(rstpm2_scenario_id = row_number()) %>%
  relocate(rstpm2_scenario_id)

scenarios_rstpm

for(i in 1:nrow(scenarios_rstpm)){
  if(!dir.exists(paste0("scen_rstpm2_", i))){
    dir.create(paste0("scen_rstpm2_", i))
  }  
}


#######################################################
# create tibble/data.frame of parameters to call from
#######################################################

set.seed(132497059)

pars_rstpm <- expand_grid(rstpm2_scenario_id = scenarios_rstpm$rstpm2_scenario_id,
                          isim = seq(nsim)) %>%
  left_join(scenarios_rstpm, by = join_by(rstpm2_scenario_id)) %>%
  mutate(data_file = paste0("../dgm", dgm, "/sim_data", isim ,".rds")) %>%
  mutate(save_file = paste0("../scen_rstpm2_", rstpm2_scenario_id, "/rstpm2_est", isim,".rds")) %>%
  mutate(seed = round(runif(n(), min = 10*7, max = 10^9))) %>%
  select(-c(rstpm2_scenario_id, dgm, isim)) 

#summary(as.factor(pars_rstpm$data_file))
#summary(as.factor(pars_rstpm$dgm))
head(pars_rstpm)

#length(unique(pars_flex$save_file))
#length(unique(pars_flex$data_file))
# test flex_fit

temp <- readRDS("scen_flex1/flex_est1.rds",)

rstpm2_fit(df_base = 5,
           df_tvc = 5, 
           data_file ="dgm2/sim_data10.rds",
           save_file = "test/rstpm_test.rds", 
           seed = 100)

test2 <- readRDS("test/rstpm_test.rds")
test
test2
#temp
#temp

rstpm2_fit <- function(df_base, df_tvc, data_file, save_file, seed){
  
  set.seed(seed) 
  
  data <- readRDS(data_file)
  
  data <- data %>%
    mutate(trt = as.numeric(trt))
  
  stpm2mod0 <- try(stpm2(Surv(time, event)~trt, 
                         data = data, 
                         df=df_base,
                         tvc = list(trt = df_tvc)))
  
  mod0_run_status <- inherits(stpm2mod0, "try-error")
  
  if(!mod0_run_status){
    
    est0 <- predict(stpm2mod0, type = "rmstdiff", 
                    newdata = tibble(trt = 0, time=5),
                    var = "trt", se.fit = TRUE )
    
    est_irmst <- tibble(estimand = "irmst", 
                        trt = NA, # needs to be 0,1.
                        t = 5,
                        value = -est0$Estimate , 
                        value_ci_low = -est0$upper,
                        value_ci_high = -est0$lower, 
                        value_se = (est0$upper-est0$lower)/(2*1.96)) %>%
      select(estimand, trt, t, value, value_ci_low, value_ci_high, value_se)
    
    est <- est_irmst
    
    saveRDS(est, save_file)
  }
}

sjob_rstpm2 <- slurm_apply(rstpm2_fit , 
                           pars_rstpm, 
                           jobname ="rstpm2_trt",
                           nodes = 10, 
                           cpus_per_node = 4, 
                           submit = TRUE,
                           pkgs = c("dplyr", "tidyr", "ggplot2", "readr",
                                    "flexsurv", "survextrap", "rstan", "emg",
                                    "pracma","rstpm2"),
                           slurm_options = list(time='02:00:00',
                                                partition='core',
                                                "mem-per-cpu"= '1G'))
system(check_status)

##########################################################


for(i in scenarios_rstpm$rstpm2_scenario_id){
  
  est <- NULL
  
  for(j in 1:nsim){
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
    if(file.exists(paste0("scen_rstpm2_", i, "/rstpm2_est", j,".rds"))){
      temp <- readRDS(paste0("scen_rstpm2_", i, "/rstpm2_est", j,".rds"))  %>%
        filter(estimand == "irmst") %>%
        mutate(isim = j, mod_type = "sim") 
      
      
      est <- rbind(est, temp)
    }
    
    
  }
  #nrow(est)
  
  saveRDS(est, paste0("scen_rstpm2_", i, "_est_rmst.rds"))
  
  true <- readRDS(paste0("dgm",scenarios_rstpm$dgm[i],"/true_est.RDS"))
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true") %>%
    filter(t == 5) %>%
    filter(estimand == "irmst")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("scen_rstpm2_", i, "_res_rmst.rds"))
  
}
