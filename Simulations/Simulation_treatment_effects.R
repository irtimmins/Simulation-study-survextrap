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

#############################################################
# specify jobname and path where results are to be stored.
#############################################################

# set to project root.
# setwd("/home/klvq491/Documents/simulation_survextrap/simsurvextrap/")

jobname <- "trt8_time" 
user <- Sys.info()["user"]
scratch_directory <- paste0("/wscratch/", user, "/")

store_res <- paste0(scratch_directory, "simsurvextrap_slurm_", jobname, "/") 

store_sjob <- paste0(store_res, "sjob_", jobname, ".rds") 
store_sjob_res <- paste0(store_res, "res_sjob_", jobname, ".csv") 
store_scenarios <- paste0(store_res,"scenarios_", jobname, ".csv") 
store_scenarios_nsim <- paste0(store_res,"scenarios_nsim_", jobname, ".csv") 
store_parameters <- paste0(store_res,"parameters_", jobname, ".csv")
store_plot_res <- paste0(store_res,"plots/","plot_res_", jobname, ".rds") 


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

if(!dir.exists(store_res)){
  dir.create(store_res)
}
setwd(store_res)


# create log.
if(!dir.exists("log")){
  dir.create("log")
}  

#savehistory(file = paste0(store_res, "/log/", jobname, "_output.txt"))

#########################################################
# DGM: import trial data
#########################################################

# from cetux control arm
surv_df <- cetux[cetux$treat=="Cetuximab",]

# end of trial follow-up time for rmst computation

maxT_data <- max(surv_df$years)

k_true <-  3
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)

#######################################################
# Specify survextrap scenarios
########################################################

# write scenarios as tibble

nsim <- 50 # simulation replicates

set.seed(378271147)

# create scenarios tibble with dgm, prior and fit.

scenarios <- expand_grid(dgm_mod = "true_mod",
                    N = 400, # trial sample size
                    maxT = maxT_data,
                    hr_scenario = 1:3) %>%
  mutate(dgm_seed = round(runif(n(), min = 10*7, max = 10^9)),
         dgm_id = row_number()) %>%
  expand_grid( prior_hscale_mean = 0,
               prior_hscale_sd = 20,
               prior_hsd_shape = 2,
               prior_hsd_rate = 1,
               prior_hrsd_shape = c(2, NA),
               prior_hrsd_rate = c(1, NA), # c(1,2,5,10,20,40, NA)
               mspline_degree = 3,
               mspline_bsmooth = c(TRUE),
               mspline_df = c(3,6,10),
               hscale_default = TRUE,
               smooth_model = c("random_walk"),
               facet_plot = TRUE,
               stan_fit_method = c("mcmc", "opt"),
               chains = c(4, NA),
               iter = c(2000, NA),
               nonprop = c(F,T)) %>% 
  filter(mspline_bsmooth == TRUE | mspline_df > 3) %>%
  filter(stan_fit_method == "opt" | !(is.na(chains) | is.na(iter))) %>% # filter out redundant stan_fit scenarios.
  filter(stan_fit_method == "mcmc" | (is.na(chains) & is.na(iter))) %>%
  filter(nonprop == T | (is.na(prior_hrsd_rate) & is.na(prior_hrsd_shape)), # filter redundant PH scenarios
         nonprop == F | !(is.na(prior_hrsd_rate) | is.na(prior_hrsd_shape))) %>%
  filter(prior_hrsd_rate == 1 | mspline_df == 10) %>%
  mutate(scenario_id = row_number()) %>%
  mutate(label_df = paste0("df~`=`~", mspline_df),
         label_gamma = paste0("sigma", "~`~`~", "`G`*`amma`*`(`*", 
                              prior_hsd_shape, "*`,`*" , 
                              prior_hsd_rate, "*`)`"),
         label_eta = paste0("eta", "~`~`~", "`N`*`(`*", 
                            prior_hscale_mean, "*`,`*" , 
                            prior_hscale_sd, "*`)`"),
         label_hrsd = paste0("`G`*`amma`*`(`*", 
                              prior_hrsd_shape, "*`,`*" , 
                              prior_hrsd_rate, "*`)`"),
         label_method = stan_fit_method) %>%
  mutate(label_hrsd = case_when(is.na(prior_hrsd_rate) == TRUE ~ "PH model",
                                is.na(prior_hrsd_rate) == FALSE ~ label_hrsd)) %>%
  mutate(label_df = fct_inorder(label_df),
         label_gamma = fct_inorder(label_gamma),
         label_eta = fct_inorder(label_eta),
         label_hrsd = fct_inorder(label_hrsd),
         label_method = fct_inorder(label_method))

View(scenarios)
#summary(as.factor(scenarios$label_hrsd))
#View(scenarios)
#scenarios <- read_csv("/wscratch/klvq491/simsurvextrap_slurm_cetux_trt1_opt/scenarios_cetux_trt1_opt.csv")
#scenarios$N <- 400
write_csv(scenarios, file = store_scenarios)
#View(scenarios)
nscen <- nrow(scenarios)

for(i in 1:nscen){

  scenario_prior <-  list(prior_hscale = p_normal(scenarios$prior_hscale_mean[i],
                                                  scenarios$prior_hscale_sd[i]),
                          prior_hsd = p_gamma(scenarios$prior_hsd_shape[i],
                                              scenarios$prior_hsd_rate[i]))
  
  if(scenarios$nonprop[i] == TRUE) {
    scenario_prior <- append(scenario_prior,
                               list(prior_hrsd = p_gamma(scenarios$prior_hrsd_shape[i],
                                                         scenarios$prior_hrsd_rate[i])))
    
  
  }
  
  scenario_mspline <- list(degree = scenarios$mspline_degree[i],
                           bsmooth = scenarios$mspline_bsmooth[i],
                           df = scenarios$mspline_df[i])
  
  scenario_stan_fit <- list(fit_method = scenarios$stan_fit_method[i])
  
  assign(paste0("prior",i), scenario_prior)
  assign(paste0("mspline",i), scenario_mspline)
  assign(paste0("stan_fit",i), scenario_stan_fit)
  

  assign(paste0("chains",i), scenarios$chains[i])
  assign(paste0("iter",i), scenarios$iter[i])
  
  assign(paste0("smooth_model", i), scenarios$smooth_model[i])
}


########################################################
# simulate datasets from dgm
########################################################

# log hazard ratio coefficients.
beta1 <- log(0.7)
beta2 <- c(log(0.5)/(1-tanh(-1.2)), 0.8, -1.2)
beta3 <- c(0, -2.8, 0.8, 0.4, 0.35)

# tibble with number of distinct data generating mechanisms.
scenarios_dgm <- scenarios %>%
  filter(!duplicated(dgm_id)) %>%
  dplyr::select(dgm_mod, hr_scenario, N, maxT, dgm_seed)

ndgm <- nrow(scenarios_dgm)

# create dgm folders.

for(i in 1:nrow(scenarios_dgm)){
  if(!dir.exists(paste0("dgm", i))){
    dir.create(paste0("dgm", i))
  }  
}

# memoise the sim_dgm_trt function so that we don't need to re-simulate data if it is already simulated under the same conditions
sim_dgm_trt_memoise <- memoise::memoise(sim_dgm_trt, cache = cachem::cache_disk(dir = paste0(store_res, ".cache")))

gen_dgm <- function(hr_function, i, N, maxT = 100, save_file, return_model=FALSE, seed ){
  
  res <- sim_dgm_trt_memoise(gamma = true_mod$coefficients,
                             knots = true_mod$knots,
                             hr_scenario = hr_function,
                             beta = get(paste0("beta", hr_function)),
                             N = N,
                             nsim = 1, 
                             maxT = maxT,
                             seed = seed,
                             nodes = 50) 
  res$i <- i # label nsim iteration number
  if(!is.null(save_file)) saveRDS(res, file = save_file)
  if(return_model == TRUE)  res
}

set.seed(378271147)

pars_dgm <- expand_grid(hr_function = 1:3,
                        i = seq(nsim)) %>%
  mutate(N = scenarios_dgm$N[1],
         maxT = scenarios_dgm$maxT[1]) %>%
  mutate(save_file = paste0("../dgm", hr_function, "/sim_data",i,".rds"),
         return_model = FALSE) %>%
  mutate(seed = sample(10^6:10^8, n(), replace= FALSE)) %>%
  slice(sample(1:n())) 


#################################################################
## Generate data from each dgm using slurm script.
#################################################################


sjob_dgm <- slurm_apply(gen_dgm, 
                        pars_dgm, 
                        jobname = "dgm",
                        nodes = 40, 
                        cpus_per_node = 4, 
                        submit = TRUE,
                        global_objects = c("true_mod",
                                           paste0("beta",1:3),
                                           "sim_dgm_trt_memoise"),
                        pkgs = c("dplyr", "tidyr", "ggplot2", "readr",
                                 "flexsurv", "survextrap", "rstan", "emg",
                                 "pracma","rstpm2" ),
                        slurm_options = list(time='02:00:00',
                                             partition='core',
                                             "mem-per-cpu"= '16G'))

check_status <- paste0("sacct -S ", as.character(Sys.Date()-20),
                       " -u ", user ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")

system(check_status)

#################################################################
## Combine results in each dgm folders into single data frame
#################################################################

for(i in 1:nrow(scenarios_dgm)){
  comb_block <- NULL
  comb <- NULL
  
  # combine datasets in blocks of 100.
  for(j in 1:max(pars_dgm$i)){
    
    block_size <- min(100, max(pars_dgm$i)) # if nsim < 100 then block_size = nsim. 
    
    if((j %% block_size) == 1) { comb <- readRDS(paste0("dgm", i, "/sim_data", j, ".rds")) 
    }    else if ((j %% block_size) != 1) { 
      temp <-  readRDS(paste0("dgm", i, "/sim_data", j, ".rds"))
      comb <- rbindlist(list(comb, temp))
    }     
    
    if((j %% block_size) == 0)  { 
    print(paste0(j, "/", nsim, ", scenario ", i)) 
    temp_block <- data.frame(comb)
    comb_block <- rbindlist(list(comb_block, temp_block))
    }
 
  }  
  
  sim_data_all <- data.frame(comb_block)
  assign(paste0("sim_data",i), sim_data_all)
  saveRDS(sim_data_all, paste0("dgm", i, "/sim_data_all.rds"))
 
}

#summary(sim_data1$i)
#summary(as.factor(sim_data1$id))
#summary(as.factor(sim_data1$i))
#max(sim_data1$i)

#test <- sim_data1 %>%
# filter(i == 10)
#nrow(test)
#head(test)
#nrow(sim_data1)
#summary(sim_data1)
#summary(as.factor(sim_data1$id))

############################################
# create directories to store results
############################################

for(i in 1:nscen){
  if(!dir.exists(paste0("scen", i))){
    dir.create(paste0("scen", i))
  }
}

#######################################################
# create tibble/data.frame of parameters to call from
#######################################################

t_vec <- seq(from = 0, to = maxT_data, length.out = 100)
rmst_vec <- 5

set.seed(132497059)

scenarios_nsim <- expand_grid(scenario_id = scenarios$scenario_id,
                              isim = 1:nsim) %>%
  left_join(scenarios %>% select(scenario_id, dgm_id, nonprop, hr_scenario), 
            by = join_by(scenario_id)) %>%
  mutate(prior_chr = paste0("prior", scenario_id),
         mspline_chr = paste0("mspline", scenario_id), 
         stan_fit_chr = paste0("stan_fit", scenario_id), 
         chains_chr = paste0("chains", scenario_id), 
         iter_chr = paste0("iter",scenario_id), 
         smooth_model_chr = paste0("smooth_model", scenario_id),
         hr_function = hr_scenario,
         maxT = maxT_data,
         sim_data_chr = paste0("sim_data",dgm_id),
         save_file = paste0("../scen", scenario_id, "/est", isim, ".rds"),
         seed =  round(runif(n(), min = 10*7, max = 10^9))) %>%
  mutate(num = row_number()) %>%
  select(-c(hr_scenario))

write_csv(scenarios_nsim, file = store_scenarios_nsim)

pars <- scenarios_nsim %>%
  select(-c(scenario_id, dgm_id, num)) %>%
  relocate(sim_data_chr) %>%
  mutate(num = row_number()) %>%
  slice(sample(1:n())) ## randomly sorting rows to even out runtime across arrays,
## increasing efficiency
#names(pars)  
write_csv(pars, file = store_parameters)
pars <- read_csv(file = store_parameters)
##############################################################
# fit survextrap models to each of the scenarios / data.frames
##############################################################

fit_est_slurm <- function(sim_data_chr,
                          isim,
                          nonprop,
                          prior_chr,
                          mspline_chr,
                          stan_fit_chr,
                          chains_chr,
                          iter_chr,
                          smooth_model_chr,
                          hr_function,
                          maxT,
                          save_file,
                          seed,
                          num){
  
  start_all <- Sys.time()
  
  set.seed(seed)
  
  sim_data <- get(sim_data_chr)
  prior <- get(prior_chr)
  mspline <- get(mspline_chr)
  stan_fit <- get(stan_fit_chr)
  smooth_model <- get(smooth_model_chr)
  print(smooth_model)
  
  if(stan_fit == "mcmc"){  
    chains <- get(chains_chr)
    iter <- get(iter_chr)
  } else {
    chains <- iter <- NA
  }
  
  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  print(nrow(data))
  
  start_fit <- Sys.time()
    
    fitted_mod <- survextrap(Surv(time, event) ~ trt,
                             nonprop = nonprop,
                             mspline = mspline, 
                             prior_hscale = prior$prior_hscale,
                             prior_hsd = prior$prior_hsd,
                             prior_hrsd = prior$prior_hrsd,
                             data = data,
                             fit_method=stan_fit$fit_method,
                             chains = chains,
                             iter = iter,
                             smooth_model = smooth_model)  


  end_fit <- Sys.time()

  print(fitted_mod)
  
  diags_df <- diags(fitted_mod) %>% 
    rename(estimand=name) 
  
  start_est <- Sys.time()
  
  est_df <- mod_est_trt(fitted_mod, estimand = NULL, t = t_vec, rmst_t = rmst_vec, nonprop = FALSE,
          return_model = TRUE, save_file = save_file)

  end_est <- Sys.time()
 
  end_all <- Sys.time()
  
  timings <- tibble("estimand" = c("run_time", "fit_time", "est_time"),
                    "value" = c(as.numeric(difftime(end_all, start_all, units='mins')),
                                as.numeric(difftime(end_fit, start_fit, units='mins')),
                                as.numeric(difftime(end_est, start_est, units='mins'))
                                )  )
    
  res <- est_df %>% 
    full_join(diags_df, by=c("estimand","value")) %>%
    full_join(timings, by=c("estimand","value") )
  
  print(tail(res)) 
  
  if (!is.null(save_file)) saveRDS(res, file=save_file)
  
  
}

## memoise fit_est_slurm
fit_est_slurm_memoise <- memoise::memoise(fit_est_slurm, 
                                          cache = cachem::cache_disk(dir = paste0(store_res, "fit.cache")))


#######################################
# # test iteration.
#######################################

#debugonce(fit_est_slurm)
#undebug(fit_est_slurm)
#i <- 150
#View(pars)
fit_est_slurm_memoise(sim_data_chr = pars$sim_data_chr[i],
              isim = pars$isim[i],
              nonprop = pars$nonprop[i],
              prior_chr = pars$prior_chr[i],
              mspline_chr = pars$mspline_chr[i],
              stan_fit_chr = pars$stan_fit_chr[i],
              chains_chr = pars$chains_chr[i],
              iter_chr = pars$iter_chr[i],
              smooth_model_chr = pars$smooth_model_chr[i],
              hr_function = pars$hr_function[i],
              maxT = pars$maxT[i],
              save_file = NULL,
              seed = pars$seed[i],
              num = pars$num[i])

###
#purrr::pmap(pars %>% slice(100), fit_est_slurm)
###

######################################
## slurm job for model fitting
######################################

# for survextrap
# nrow(pars)
# nsim = 50 ->, 100 arrays takes < 30 mins each.
# nsim = 1,000 -> 100 arrays takes ~12 hours each

attach_obj <- c("surv_df", 
                paste0("sim_data",1:ndgm),
                paste0("prior",1:nscen),
                paste0("mspline",1:nscen),
                paste0("stan_fit",1:nscen),
                paste0("chains",c(1:nscen)),
                paste0("iter",c(1:nscen)),
                paste0("smooth_model",1:nscen),
                "maxT_data",
                "fit_model", "fit_model_slurm",
                "catchToList", "throw_sampler_warnings",
                "extract_data",
                "rmst_vec", "t_vec",
                "mod_est_trt", "true_mod", "diags")

attach_pkgs <-  c("dplyr", "tidyr", "ggplot2", "readr",
                  "flexsurv", "survextrap", "rstan")

slurm_opts <- list(time='23:00:00',
                   partition='core',
                   "mem-per-cpu"= '16G')

# nsim=50 with N=400 takes ~30 mins.
# nsim=10 with N=10,000 opt takes ~10 mins.

sjob_fit_est <- slurm_apply(fit_est_slurm_memoise, 
                            pars, 
                            jobname = "fit_est",
                            nodes = 80, 
                            cpus_per_node = 4, 
                            submit = TRUE,
                            global_objects = attach_obj,
                            pkgs = attach_pkgs,
                            slurm_options = slurm_opts)

saveRDS(sjob_fit_est, file = store_sjob)
system(check_status)


test <- readRDS("scen3/est6.rds")
summary(as.factor(test$estimand))
nrow(test)
View(test)


##########################################################
# Get true HR/RMST etc. for each dgm.
##########################################################

# Generate big datasets.
set.seed(378271147)
pars_dgm_big <- expand_grid(hr_function = 1:3,
                        i = 1:200) %>%
  mutate(N = 1e5,
         maxT = 100) %>%
  mutate(save_file = paste0("../dgm", hr_function, "/sim_data_big",i,".rds"),
         return_model = FALSE) %>%
  mutate(seed = sample(10^6:10^8, n(), replace= FALSE)) %>%
  slice(sample(1:n())) # randomise rows to even out slurm run time across arrays.

# sum(pars_dgm$N)
# length(unique(pars_dgm$seed))
# each array takes ~1 hour.

sjob_dgm_big <- slurm_apply(gen_dgm, 
                        pars_dgm_big, 
                        jobname ="dgm_big",
                        nodes = 60, 
                        cpus_per_node = 4, 
                        submit = TRUE,
                        global_objects = c("true_mod",
                                           paste0("beta",1:3),
                                           "sim_dgm_trt_memoise"),
                        pkgs = c("dplyr", "tidyr", "ggplot2", "readr",
                                 "flexsurv", "survextrap", "rstan", "emg",
                                 "pracma","rstpm2" ),
                        slurm_options = list(time='02:00:00',
                                             partition='core',
                                             "mem-per-cpu"= '16G'))

system(check_status)

## Combine results into big data frame with N = 10^7 patients in each arm.

for(i in 1:nrow(scenarios_dgm)){
  for(j in 1:max(pars_dgm_big$i)){
    if(j == 1) { comb <- readRDS(paste0("dgm", i, "/sim_data_big", j, ".rds")) }
    else { 
      temp <-  readRDS(paste0("dgm", i, "/sim_data_big", j, ".rds"))
      comb <- rbindlist(list(comb, temp))
    }
    if((j %% 10) == 0) print(paste0(j, "/",max(pars_dgm_big$i), " dgm ", i)) # track progress.
  }  
  sim_data_big <- data.frame(comb)
  saveRDS(sim_data_big, paste0("dgm", i, "/sim_data_big.rds"))
}


####################################################
# Evaluate true values for irmst and hazard ratio. 
####################################################

# irmst is the incremental rmst (gain in rmst), 
# taken across contrasts trt = 0/1.

for(i in 1:nrow(scenarios_dgm)){
  
  hr_scenario <- scenarios$scenario_id[i]
  
  # get hr and haz functions called from "functions_treatment_effect.R".
  
  hr <- scenario_hr(hr_scenario)
  haz <- scenario_haz(hr_scenario)
  
  beta_true <- get(paste0("beta", scenarios$scenario_id[i]))
  
  # haz_true contains the true hazard rate
  # for each time point in t_vec.
  # this is evaluated analytically, since the hazard 
  # function for each arm has a closed-form.
  
  haz_true <- expand_grid(estimand = "hazard",
                          trt = 0:1,
                          t = t_vec,
                          value = NA,
                          value_ci_low = NA,
                          value_ci_high = NA,
                          value_se = NA ) %>%
    arrange(t)
  
  for(j in 1:length(t_vec)){
    
    haz_true$value[haz_true$t == t_vec[j]] <- haz(t =  t_vec[j],
                                                  gamma = true_mod$coefficients,
                                                  knots = true_mod$knots, 
                                                  beta =  beta_true, 
                                                  trt = 0:1) 
    
    
  }
  
  # hr_true contains the true hazard ratio 
  # each time point in t_vec.
  # this is evaluated analytically, since hr
  # has a closed-form.
  
  hr_true <- tibble(estimand = "hr",
                    trt = NA,
                    t = t_vec) %>%
    mutate(value = hr(t = t_vec,
                     gamma = true_mod$coefficients,
                     knots = true_mod$knots, 
                     beta =  beta_true) ) %>%
    mutate(value_ci_low = NA,
           value_ci_high = NA,
           value_se = NA)

  sim_data_big <- readRDS(paste0("dgm", i, "/sim_data_big.rds"))
  
  # rmst_true and irmst_true contains the true rmst and irmst (incremental rmst)
  # for each time point in rmst_vec.
  # As this is not analytically tractable,
  # the truth' is estimated from the large data set 'sim_data_big'
  # that was created above.
  
  rmst_true <- expand_grid(estimand = "rmst",
                           trt = c(0,1),
                           t = rmst_vec,
                           value = NA,
                           value_ci_low = NA,
                           value_ci_high = NA,
                           value_se = NA )
  
  irmst_true <- tibble(estimand = "irmst",
                       trt = NA,
                       t = rmst_vec,
                       value = NA,
                       value_ci_low = NA,
                       value_ci_high = NA,
                       value_se = NA )
  
  # for each time point in rsmt_vec
  # the true irmst is computed and put in the "value" column
  
  print("rmst and irmst..")
  
  for(j in 1:length(rmst_vec)){
    
    print(paste0(j, "/",length(rmst_vec), " rmst")) # track progress.
    #  j <- 6
    rmst_temp <- sim_data_big %>%
      mutate(time2 = pmin(rmst_vec[j], time)) %>%
      group_by(trt) %>%
      summarise(rmst = mean(time2))  %>%
      mutate(irmst = rmst - lag(rmst, default = 0)) %>%
      pivot_longer(cols = c("rmst", "irmst"), names_to = "estimand", values_to = "value") %>%
      filter(trt == 1 | estimand == "rmst") %>%
      mutate(trt = case_when(
        estimand == "irmst"  ~ NA,
        .default = trt) ) 
    
    irmst_true$value[j] <- rmst_temp %>%
      filter(estimand == "irmst") %>%
      pull(value)
    
    rmst_true$value[rmst_true$t == rmst_vec[j]] <- rmst_temp %>%
      filter(estimand == "rmst") %>%
      pull(value)
      
  }
  
  # surv_true contains the true survival probabilities
  # for each time point in t_vec.
  # As this is not analytically tractable,
  # the survival probabilities are estimated from the
  # large data set 'sim_data_big' simulated above.
  
  surv_true <- expand_grid(estimand = "survival",
                           trt = 0:1,
                           t = t_vec,
                           value = NA,
                           value_ci_low = NA,
                           value_ci_high = NA,
                           value_se = NA ) %>%
    arrange(t)
  
  # for each time point in t_vec
  # the survival probability is computed and put
  # in the "value" column
  print("survival probs..")
  for(j in 1:length(t_vec)){
    
  if((j %% 10) == 0)  print(paste0(j, "/",length(t_vec), " survival probs")) # track progress.
    
  surv_true$value[surv_true$t == t_vec[j]] <-     
      sim_data_big %>%
      mutate(time2 = time > t_vec[j]) %>%
      group_by(trt) %>%
      summarise(survival = sum(time2)/n()) %>%
      select(survival) %>%
      pull()
  
    }
  
  true_est <- rbind(rmst_true,
                    irmst_true,
                    hr_true,
                    haz_true,
                    surv_true)
  
  saveRDS(true_est, file = paste0("dgm", i, "/true_est.RDS"))

}

test <- readRDS( paste0("dgm", 3, "/true_est.RDS"))
summary(as.factor(test$estimand))
View(test)

####################################################
# combine results for estimands into tidy format.
####################################################

# try without using slurm.

for(i in 1:nscen){
  
  dgm_num <- scenarios$dgm_id[i]
  est <- NULL

  for(j in 1:nsim){
    # print progress track every 10th iteration
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
    if(file.exists(paste0("scen", i, "/est", j,".rds"))){
      temp <- readRDS(paste0("scen", i, "/est", j,".rds"))
      temp <- temp %>%
        mutate(isim = j, mod_type = "sim") 
      est <-  rbindlist(list(est, temp))
    }
  }
  est <- data.frame(est)
  saveRDS(est, paste0("scen", i, "_est.rds"))
  
  # read in true values of estimands.
  true <- readRDS(paste0("dgm", dgm_num, "/true_est.RDS"))
  true <- true %>%
    mutate(isim = 0, mod_type = "true")
  
  # append true values alongside estimates.
  res <- rbind(true, est) %>%
    group_by(estimand, isim) %>%
    mutate(id = row_number(),
           estimand_id = paste0(estimand, id)) %>%
    relocate(estimand_id) %>%
    select(-id) %>%
    ungroup()
  
  saveRDS(res, paste0("scen", i, "_res.rds"))
  
}


#################################
## Check results.
#################################

i <- 2
test <- readRDS(paste0("scen", i, "_res.rds"))
summary(as.factor(res$estimand))
summary(as.factor(res$isim))

####################################################
# Create directory for storing plots.
####################################################

if(!dir.exists("plots")){
  dir.create("plots")
}  

####################################################
# Prep data for plotting hazard ratio.
####################################################

# scenarios

for(i in scenarios$scenario_id){
  
  temp <- readRDS(paste0("scen", i, "_res.rds"))
  
  temp <- temp %>% 
    mutate(scenario_id = i) # %>%

  
  if(i == scenarios$scenario_id[1]){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }
  print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
  res <- data.frame(res)
}

#scenarios <- read_csv(file = store_scenarios)
res <- res %>%
  left_join(scenarios %>% 
              select(hr_scenario, prior_hscale_mean, prior_hscale_sd,
                     prior_hsd_shape, prior_hsd_rate, 
                     prior_hrsd_shape, prior_hrsd_rate, mspline_degree, mspline_df,
                     mspline_bsmooth,smooth_model,
                     stan_fit_method, nonprop, prior_hrsd_rate, scenario_id,
                     label_df,  label_gamma, label_eta,  label_hrsd, label_method), 
            by = join_by(scenario_id))

saveRDS(res, file = store_plot_res)
#names(res)
#head(res)
res <- readRDS(store_plot_res)



