
library(tidyr)
library(dplyr)
library(ggplot2)
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
library(rstpm2)

#############################################################
# specify jobname and path where results are to be stored.
#############################################################

setwd(here::here())

jobname <- "cetux16"
jobname <- "niv8" 
user <- Sys.info()["user"]
scratch_directory <- paste0("/wscratch/", user, "/")

store_res <- paste0(scratch_directory, "simsurvextrap_slurm_", jobname, "/") 

#############################################################
# import trial data and call scripts to create functions.
#############################################################
setwd(here::here())
source("R/simulate_dgm.R")
source("R/fit_model.R")
source("R/estimands.R")

#dir.create(store_res)
setwd(store_res)

#########################################################
# DGM: import trial data
#########################################################

# From Cetux control arm
surv_df <- cetux[cetux$treat=="Cetuximab",]

# end of trial follow-up time for rmst computation

maxT_data <- max(surv_df$years)

k_true <-  3
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)

# Nivolumav pfs

surv_df <- read_csv(here::here("nivo_data/nivolumab_digiKM_pfs.csv"),
                    show_col_types = FALSE)  %>%
  filter(treat==2) %>%
  mutate(years = time/12) %>%
  rename(d = status)

# end of trial follow-up time for rmst computation

maxT_data <- max(surv_df$years)

k_true <-  6
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)

###################################################

true <- mod_est(mod = true_mod, estimand = NULL, t = 5, rmst_t = 5, nonprop = FALSE,
                return_model = TRUE, save_file = NULL)

saveRDS(true, "true_est.RDS")

set.seed(242642938)

#######################################################
# Specify survextrap scenarios
########################################################

# write scenarios as tibble

nsim <- 1000 # simulation replicates

scenarios_dgm <- tibble("dgm_mod" = "true_mod",
                        "N" = c(200), # trial sample size
                        "maxT" = maxT_data,
                        "dgm_seed" = 126283921) %>%
  mutate(dgm_id = row_number())


########################################################
# simulate dataset from dgm
########################################################

# may take a few minutes if nsim large.
ndgm <- nrow(scenarios_dgm)

for(i in 1:ndgm){
  
  true_mod <- get(scenarios_dgm$dgm_mod[i])
  
  assign(paste0("sim_data",i), sim_dgm(model_class = "flexsurv",
                                       true_mod = true_mod,
                                       N = scenarios_dgm$N[i],
                                       nsim = nsim, 
                                       maxT = scenarios_dgm$maxT[i],
                                       seed = scenarios_dgm$dgm_seed[i]) )
  
}


############################################
# create directories to store results
############################################

#getwd()

nsim <- 1000

scenarios_flex <- expand_grid(k_knots = c(0,1,2,4,8)) %>%
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

pars_flex <- expand_grid(sim_data_chr = "sim_data1",
                         flex_scenario_id = scenarios_flex$flex_scenario_id,
                         isim = seq(nsim)) %>%
  left_join(scenarios_flex, by = join_by(flex_scenario_id)) %>%
  mutate(save_file = paste0("../scen_flex", flex_scenario_id, "/flex_est", isim,".rds")) %>%
  mutate(seed =  round(runif(n(), min = 10*7, max = 10^9))) %>%
  select(-c(flex_scenario_id))  %>%
  slice(sample(1:n())) 

saveRDS(pars_flex, file = paste0("pars_flex_",jobname, ".rds"))
pars_flex <- readRDS(file = paste0("pars_flex_",jobname, ".rds"))
summary(as.factor(pars_flex$k_knots))
head(pars_flex)
#length(unique(pars_flex$save_file))
#length(unique(pars_flex$data_file))

# test flex_fit

flex_fit(sim_data_chr = "sim_data1",
         isim = 10,
         k_knots = 6,
         save_file = "scen_flex1/flex_test1.rds",
         seed = 70494313)

temp <- readRDS("scen_flex1/flex_test1.rds",)
temp
#temp

flex_fit <- function(sim_data_chr, isim, k_knots, save_file, seed){
  
  set.seed(seed) 
  
  sim_data <- get(sim_data_chr)
  
  
  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  
  flexmod0 <- try(flexsurvspline(Surv(time, event)~1, 
                                 data = data, 
                                 k = k_knots,
                                 spline = "splines2ns"))
  
  
  mod0_run_status <- inherits(flexmod0, "try-error")
  
  # model run succesfully, 
  # continue with standsurv:
  
  if(!mod0_run_status){
    
    covar_status <- !is.na(flexmod0$cov[1])
    
    if(covar_status){
      
      est0 <- standsurv(flexmod0, type = "rmst", t = 5, ci = TRUE, se = TRUE, boot = T, B = 1000)
      
      est_rmst <- est0 %>%
        mutate(estimand = "rmst", 
               trt = NA, # needs to be 0,1.
               t = 5,
               value = at1, 
               value_ci_low = at1_lci,
               value_ci_high = at1_uci, 
               value_se = at1_se) %>%
        select(estimand, trt, t, value, value_ci_low, value_ci_high, value_se)
      
      
      est <- est_rmst
      
      saveRDS(est, save_file)
    }
    
  }
}

sjob_flex <- slurm_apply(flex_fit , 
                         pars_flex, 
                         jobname ="flex",
                         nodes = 20, 
                         cpus_per_node = 4, 
                         submit = TRUE,
                         global_objects = c("sim_data1"), 
                         pkgs = c("dplyr", "tidyr","readr",
                                  "flexsurv", "survextrap", "rstan","rstpm2" ),
                         slurm_options = list(time='02:00:00',
                                              partition='core',
                                              "mem-per-cpu"= '4G'))

check_status <- paste0("sacct -S ", as.character(Sys.Date()-4),
                       " -u ", user ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")
system(check_status)
temp <- readRDS("../simsurvextrap_slurm_master_cetux12/_rslurm_flex/results_19.RDS")
#temp <- readRDS("../simsurvextrap_slurm_cetux_trt4/_rslurm_flex/results_14.RDS")
temp


############################################
# create directories to store results
############################################

#getwd()

nsim <- 1000

scenarios_rstpm <- expand_grid(df = c(1,2,3,5,9)) %>%
  mutate(rstpm2_scenario_id = row_number()) %>%
  relocate(rstpm2_scenario_id)

for(i in 1:nrow(scenarios_rstpm)){
  if(!dir.exists(paste0("scen_rstpm2_", i))){
    dir.create(paste0("scen_rstpm2_", i))
  }  
}

rstpm2_fit(sim_data_chr = "sim_data1",
           isim = 10,
           df = 6,
           save_file = "scen_rstpm2_1/rstpm2_test1.rds",
           seed = 70494313)


rstpm2_fit <- function(sim_data_chr, isim, df, save_file, seed){
  
  set.seed(seed) 
  
  sim_data <- get(sim_data_chr)
  
  
  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  
  stpm2mod0 <- try(stpm2(Surv(time, event)~1, 
                         data = data , 
                         df=df))
  
  
  mod0_run_status <- inherits(stpm2mod0, "try-error")
  
  # model run succesfully, 
  # continue with standsurv:
  
  if(!mod0_run_status){
    
    est0 <- predict(stpm2mod0, type = "rmst", newdata = tibble(time=5),se.fit = TRUE)
    
    est_rmst <- tibble(estimand = "rmst", 
                       trt = NA, # needs to be 0,1.
                       t = 5,
                       value = est0$Estimate , 
                       value_ci_low = est0$lower,
                       value_ci_high = est0$upper, 
                       value_se = NA) %>%
      select(estimand, trt, t, value, value_ci_low, value_ci_high, value_se)
    
    est <- est_rmst
    print(head(est))
    saveRDS(est, save_file)
  }
}

pars_rstpm2 <- expand_grid(sim_data_chr = "sim_data1",
                           rstpm2_scenario_id =  scenarios_rstpm$rstpm2_scenario_id,
                           isim = seq(nsim)) %>%
  left_join(scenarios_rstpm, by = join_by(rstpm2_scenario_id)) %>%
  mutate(save_file = paste0("../scen_rstpm2_", rstpm2_scenario_id, "/rstpm2_est", isim,".rds")) %>%
  mutate(seed =  round(runif(n(), min = 10*7, max = 10^9))) %>%
  select(-c(rstpm2_scenario_id)) 

saveRDS(pars_rstpm2, file = paste0("pars_rsptm2_",jobname, ".rds"))
pars_rstpm2 <- readRDS(file = paste0("pars_rsptm2_",jobname, ".rds"))
#summary(as.factor(pars_rstpm2$df))
#head(pars_flex)

sjob_rstpm2 <- slurm_apply(rstpm2_fit , 
                           pars_rstpm2, 
                           jobname ="rstpm2",
                           nodes = 10, 
                           cpus_per_node = 4, 
                           submit = TRUE,
                           global_objects = c("sim_data1"), 
                           pkgs = c("dplyr", "tidyr", "ggplot2", "readr",
                                    "flexsurv", "survextrap", "rstan", "emg",
                                    "pracma","rstpm2" ),
                           slurm_options = list(time='01:00:00',
                                                partition='core',
                                                "mem-per-cpu"= '4G'))
system(check_status)

