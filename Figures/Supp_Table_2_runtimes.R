###########
library(here)
library(tidyr)
library(dplyr)
library(ggplot2)
library(survextrap)
library(survminer)
library(flexsurv)
library(rslurm)
library(rsimsum)
library(stringr)
library(readr)
library(rstan)
library(tibble)
library(purrr)
library(forcats)
library(cowplot)
library(gt)
library(kableExtra)
library(formattable)
library(data.table)

user <- Sys.info()["user"]
store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")


cetux_jobname <-   "cetux17_time"
niv_jobname <- "niv9_time" 
#jobname <- cetux_jobname

for(jobname in c(cetux_jobname, niv_jobname)){
  #jobname <-  "niv9_time" 
  #niv_jobname <- "niv9_time" 
  #jobname <- "master_niv3" 
  
  setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))
  print(getwd())
  #dir.create("plots")
  
  scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
  nscen <- nrow(scenarios)
  
  #res_test <- readRDS("scen66_est_rmst.rds")
  #res_test
  #View(res_test)
  #summary(as.factor(res_test$estimand_id))
  
  ###################################################
  # Combine data together.
  ###################################################
  
  for(i in scenarios$scenario_id){
    
   for(j in 1:100){ 
     #i <- j <- 1
     temp <- readRDS(paste0("scen", i ,"/", "est",j, ".rds"))
     if(j == 1){
      res <- temp
     }  else  {
       res <- rbindlist(list(res, temp))
     }
   }
    saveRDS(res, paste0("scen", i ,"/est_all.rds"))
  }
  
  
  for(i in scenarios$scenario_id){
    #i <- 1
    temp <- readRDS(paste0("scen", i ,"/est_all.rds"))
    
    temp <- temp %>% 
      mutate(scenario_id = i)  %>%
      filter(estimand == "fit_time")
    
    if(i == scenarios$scenario_id[1]){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
    res <- data.frame(res)

  }
  
  #names(scenarios)
  #scenarios <- read_csv(file = store_scenarios)
  fit_all <- res %>%
    group_by(scenario_id) %>%
    summarise(mean = 60*mean(value)) %>%
    left_join(scenarios %>%
                select(scenario_id, stan_fit_method, prior_hsd_rate, 
                       mspline_df, mspline_bsmooth, smooth_model), 
              join_by(scenario_id)) %>%
    filter(prior_hsd_rate == 1,
           mspline_bsmooth == T,
           smooth_model == "random_walk")
  
 #View(fit_all)
  
 #cetux_test <- fit_all
 #niv_test <- fit_all
 
  saveRDS(fit_all, "fit_time.rds")
}



fit_cetux <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", cetux_jobname, "/fit_time.rds"))
fit_niv <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", niv_jobname, "/fit_time.rds"))

fit_cetux
fit_niv


table_df <- rbind(fit_cetux %>% mutate(case_study = "Cetuximab OS"), 
                  fit_niv%>% mutate(case_study = "Nivoumab PFS")) %>%
  relocate(case_study) %>%
  select(-c(scenario_id, prior_hsd_rate, mspline_bsmooth, smooth_model)) %>%
  rename(time = mean)


table_df$time[table_df$stan_fit_method == "mcmc"] <-   sprintf("%.1f", round(table_df$time[table_df$stan_fit_method == "mcmc"] , 2))
table_df$time <- as.numeric(table_df$time)
table_df$time[table_df$stan_fit_method == "opt"] <-   sprintf("%.2f", round(table_df$time[table_df$stan_fit_method == "opt"] , 2))
table_df


write_delim(table_df, paste0(store_directory, "tables/single_arm_runtime.csv"), delim = ",")


#######################################
# Treatment effect run time.
#######################################


trt_jobname <-   "trt8_time"

#jobname <- cetux_jobname


#jobname <-  "niv9_time" 
#niv_jobname <- "niv9_time" 
#jobname <- "master_niv3" 

setwd(paste0(store_directory , "simsurvextrap_slurm_", trt_jobname, "/"))
print(getwd())
#dir.create("plots")

scenarios <- read_csv(paste0("scenarios_", trt_jobname ,".csv"))
nscen <- nrow(scenarios)

#res_test <- readRDS("scen66_est_rmst.rds")
#res_test
#View(res_test)
#summary(as.factor(res_test$estimand_id))

###################################################
# Combine data together.
###################################################

for(i in scenarios$scenario_id){
  
  for(j in 1:50){ 
    #i <- j <- 1
    temp <- readRDS(paste0("scen", i ,"/", "est",j, ".rds"))
    if(j == 1){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
  }
  saveRDS(res, paste0("scen", i ,"/est_all.rds"))
}


for(i in scenarios$scenario_id){
  #i <- 1
  temp <- readRDS(paste0("scen", i ,"/est_all.rds"))
  
  temp <- temp %>% 
    mutate(scenario_id = i)  %>%
    filter(estimand == "fit_time")
  
  if(i == scenarios$scenario_id[1]){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }
  print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
  res <- data.frame(res)
  
}

#names(scenarios)
#scenarios <- read_csv(file = store_scenarios)
fit_all <- res %>%
  group_by(scenario_id) %>%
  summarise(mean = 60*mean(value)) %>%
  left_join(scenarios %>%
              select(scenario_id,hr_scenario, stan_fit_method, nonprop, prior_hsd_rate, 
                     mspline_df, mspline_bsmooth, smooth_model), 
            join_by(scenario_id)) %>%
  filter(prior_hsd_rate == 1,
         mspline_bsmooth == T,
         smooth_model == "random_walk") %>%
  select(-c(prior_hsd_rate, mspline_bsmooth, smooth_model, scenario_id)) %>%
  relocate(hr_scenario) %>%
  arrange(nonprop) %>%
  arrange(stan_fit_method) %>%
  arrange(hr_scenario) %>%
  rename(time = mean)
 
#View(fit_all)
#View(scenarios)
#names(scenarios)
#cetux_test <- fit_all
#niv_test <- fit_all

saveRDS(fit_all, "fit_time.rds")

table_trt_df <- fit_all

table_trt_df$time[table_trt_df$stan_fit_method == "mcmc"] <-   sprintf("%.1f", round(table_trt_df$time[table_trt_df$stan_fit_method == "mcmc"] , 2))
table_trt_df$time <- as.numeric(table_trt_df$time)
table_trt_df$time[table_trt_df$stan_fit_method == "opt"] <-   sprintf("%.2f", round(table_trt_df$time[table_trt_df$stan_fit_method == "opt"] , 2))
table_trt_df

write_delim(table_trt_df, paste0(store_directory, "tables/two_arm_runtime.csv"), delim = ",")

#############################

