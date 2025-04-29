
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

#############################################################
# specify jobname and path where results are to be stored.
#############################################################

setwd(here::here())

#jobname <- "cetux17_time"
jobname <- "niv9_time" 

user <- Sys.info()["user"]
scratch_directory <- paste0("/wscratch/", user, "/")

store_res <- paste0(scratch_directory, "simsurvextrap_slurm_", jobname, "/") 

store_sjob <- paste0(store_res, "sjob_", jobname, ".rds")
store_sjob_res <- paste0(store_res, "res_sjob_", jobname, ".csv") 

store_scenarios <- paste0(store_res,"scenarios_", jobname, ".csv") 
store_scenarios_nsim <- paste0(store_res,"scenarios_nsim_", jobname, ".csv") 
store_parameters <- paste0(store_res,"parameters_", jobname, ".csv")

#############################################################
# import trial data and call scripts to create functions.
#############################################################

setwd(here::here())
source("R/simulate_dgm.R")
source("R/fit_model.R")
source("R/estimands.R")

dir.create(store_res)
setwd(store_res)

#########################################################
# DGM: import trial data
#########################################################

# Cetuximab OS

surv_df <- cetux[cetux$treat=="Cetuximab",]

# end of trial follow-up time for rmst computation

maxT_data <- max(surv_df$years)

k_true <-  3
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)


# Nivolumav PFS

surv_df <- read_csv(here::here("nivo_data/nivolumab_digiKM_pfs.csv"),
                    show_col_types = FALSE)  %>%
  filter(treat==2) %>%
  mutate(years = time/12) %>%
  rename(d = status)

# end of trial follow-up time for rmst computation

maxT_data <- max(surv_df$years)

k_true <-  6
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)


#######################################################
# Specify survextrap scenarios
########################################################

# write scenarios as tibble

nsim <- 100 # simulation replicates

scenarios_dgm <- tibble("dgm_mod" = "true_mod",
                        "N" = c(200), # trial sample size
                        "maxT" = maxT_data,
                        "dgm_seed" = 126283921) %>%
  mutate(dgm_id = row_number())

scenarios_priors <- 
  expand_grid(prior_hscale_mean = 0,
              prior_hscale_sd = 20,
              prior_hsd_shape = 2,
              prior_hsd_rate = c(1,5,20),
              mspline_degree = 3,
              mspline_bsmooth = c(TRUE,FALSE),
              mspline_df = c(3,6,10),
              hscale_default = TRUE,
              smooth_model = c("exchangeable","random_walk"),
              facet_plot = TRUE) %>%
  filter(mspline_bsmooth == TRUE | mspline_df > 3)

scenarios_fit <- tibble("stan_fit_method" = c("mcmc","opt"),
                        "chains" = c(4,NA),
                        "iter" = c(2000,NA)) %>%
  mutate(fit_id = n())

# create scenarios tibble with dgm, prior and fit.

# replicate scenarios across each stan_fit choice
scenarios <- cbind(scenarios_fit %>%
                     slice(rep(1:n(), 
                               each = nrow(scenarios_priors))),
                   scenarios_priors)
# replicate scenarios across each dgm
scenarios <- cbind(scenarios_dgm %>%
                     slice(rep(1:n(), each = nrow(scenarios))),
                   scenarios) %>%
  mutate(scenario_id = row_number())

write_csv(scenarios, file = store_scenarios)
scenarios <- read_csv(store_scenarios)
#View(scenarios)

nscen <- nrow(scenarios)

for(i in 1:nscen){
  
  scenario_prior <-  list(prior_hscale = p_normal(scenarios$prior_hscale_mean[i],
                                                  scenarios$prior_hscale_sd[i]),
                          prior_hsd = p_gamma(scenarios$prior_hsd_shape[i],
                                              scenarios$prior_hsd_rate[i]))
  
  scenario_mspline <- list(degree = scenarios$mspline_degree[i],
                           bsmooth = scenarios$mspline_bsmooth[i],
                           df = scenarios$mspline_df[i])
  
  scenario_stan_fit <- list(fit_method = scenarios$stan_fit_method[i])
  
  assign(paste0("smooth_model", i), scenarios$smooth_model[i])
  assign(paste0("prior",i), scenario_prior)
  assign(paste0("mspline",i), scenario_mspline)
  assign(paste0("stan_fit",i), scenario_stan_fit)
  
  
  assign(paste0("chains",i), scenarios$chains[i])
  assign(paste0("iter",i), scenarios$iter[i])
  
}


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

for(i in 1:nscen){
  if(!dir.exists(paste0("scen", i))){
    dir.create(paste0("scen", i))
  }
}

#######################################################
# create tibble/data.frame of parameters to call from
#######################################################

# time points for hazard and survival estimates
t_vec <- seq(from = 0, to = maxT_data, length.out = 100)
# time points for rmst estimates
rmst_vec <- 5

set.seed(132497059)

scenarios_nsim <- expand_grid(scenario_id = scenarios$scenario_id,
                              isim = 1:nsim) %>%
  left_join(scenarios %>% select(scenario_id, dgm_id), 
            by = join_by(scenario_id)) %>%
  mutate(prior_chr = paste0("prior", scenario_id),
         mspline_chr = paste0("mspline", scenario_id), 
         stan_fit_chr = paste0("stan_fit", scenario_id), 
         chains_chr = paste0("chains", scenario_id), 
         iter_chr = paste0("iter",scenario_id), 
         smooth_model_chr = paste0("smooth_model", scenario_id),
         maxT = maxT_data,
         sim_data_chr = paste0("sim_data",dgm_id),
         save_file = paste0("../scen", scenario_id, "/est", isim, ".rds"),
         seed =  round(runif(n(), min = 10*7, max = 10^9))) %>%
  mutate(num = row_number()) %>%
  slice(sample(1:n())) # randomisation of rows evens out slurm runtime.

write_csv(scenarios_nsim, file = store_scenarios_nsim)
scenarios_nsim <- read_csv(store_scenarios_nsim)

pars <- scenarios_nsim %>%
  select(-c(scenario_id, dgm_id, num)) %>%
  relocate(sim_data_chr) %>%
  mutate(num = row_number()) %>%
  slice(sample(1:n())) 

write_csv(pars, file = store_parameters)
pars <- read_csv(store_parameters)
#summary(pars)

fit_est_slurm <- function(sim_data_chr,
                          isim,
                          prior_chr,
                          mspline_chr,
                          stan_fit_chr,
                          chains_chr,
                          iter_chr,
                          smooth_model_chr,
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

  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  start_fit <- Sys.time()
  
  
  if(stan_fit == "mcmc"){  
    chains <- get(chains_chr)
    iter <- get(iter_chr)
  } else {
    chains <- iter <- NA
  }
  
  start_fit <- Sys.time()
  
  fitted_mod <- survextrap(Surv(time, event) ~ 1, 
                           mspline = mspline, 
                           prior_hscale = prior$prior_hscale,
                           prior_hsd = prior$prior_hsd,
                           data = data,
                           fit_method=stan_fit$fit_method,
                           chains = chains,
                           iter = iter,
                           smooth_model = smooth_model) 
  
  end_fit <- Sys.time()

  diags_df <- diags(fitted_mod) %>% 
    rename(estimand=name) 
  
  start_est <- Sys.time()
  
  est_df <- mod_est(fitted_mod, estimand = NULL, t = t_vec, rmst_t = rmst_vec, nonprop = FALSE,
          return_model = TRUE, save_file = NULL)
  
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
  
 #  for testing
 #  print(head(res)) 
 #  print(tail(res)) 
  
  if (!is.null(save_file)) saveRDS(res, file=save_file)
  
}

#######################################
## test iteration.
#######################################

#i <- 500
fit_est_slurm(sim_data_chr = pars$sim_data_chr[i],
              isim = pars$isim[i],
              prior_chr = pars$prior_chr[i],
              mspline_chr = pars$mspline_chr[i],
              stan_fit_chr = pars$stan_fit_chr[i],
              chains_chr = pars$chains_chr[i],
              iter_chr = pars$iter_chr[i],
              smooth_model_chr = pars$smooth_model_chr[i],
              maxT = pars$maxT[i],
              save_file = NULL,
              seed = pars$seed[i],
              num = pars$num[i])

######################################
## slurm job for model fitting
######################################

pars_slurm <- pars %>%
  separate(prior_chr, into = c("temp", "scenario_id"), sep = "prior", remove = FALSE) %>%
  select(-temp) %>% 
  relocate(scenario_id, .after = last_col()) %>%
  mutate(scenario_id = as.numeric(scenario_id)) %>%
  left_join(scenarios %>% select(scenario_id, stan_fit_method),
            by = join_by(scenario_id))  %>% 
  select(-c(stan_fit_method, scenario_id))

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
                "extract_data", "diags", 
                "rmst_vec", "t_vec",
                "mod_est", "true_mod")

attach_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr",
                 "flexsurv", "survextrap", "rstan")

slurm_opts <- list(time='23:00:00',
                   partition='core',
                   "mem-per-cpu"= '16G')

sjob_fit_est <- slurm_apply(fit_est_slurm, 
                                pars_slurm, 
                                jobname = jobname,
                                nodes = 20, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                global_objects = attach_obj,
                                pkgs = attach_pkgs ,
                                slurm_options = slurm_opts )


saveRDS(sjob_fit_est, file = store_sjob)
sjob_fit_est <- readRDS(file = store_sjob)

check_status <- paste0("sacct -S ", as.character(Sys.Date()-5),
                       " -u ", user ,
                       " --format=JobID,Jobname,partition,state,elapsed,ncpus -X")

system(check_status)

#######################################
# check slurm results.
#######################################

i_test <- 5

res <- readRDS(paste0(store_res, "_rslurm_",jobname,"/results_", i_test, ".RDS"))
res <- res[[1]]
res

####################################################
# Prepare look-up table for scenario labels.
####################################################

dgm_num <- 1

scenarios_plot <- scenarios %>%
  filter(dgm_id == dgm_num) %>%
  mutate(stan_fit_method = if_else(stan_fit_method == "mcmc", "MCMC", "Opt"))


label_lookup <- scenarios_plot %>%
  select(c(scenario_id, mspline_df, prior_hsd_shape, 
           prior_hsd_rate,prior_hscale_mean, prior_hscale_sd ,stan_fit_method)) %>%
  mutate(label_df = paste0("df~`=`~", mspline_df),
         label_gamma = paste0("sigma", "~`~`~", "`G`*`amma`*`(`*", 
                              prior_hsd_shape, "*`,`*" , 
                              prior_hsd_rate, "*`)`"),
         label_eta = paste0("eta", "~`~`~", "`N`*`(`*", 
                            prior_hscale_mean, "*`,`*" , 
                            prior_hscale_sd, "*`)`"),
         label_method = stan_fit_method) %>%
  mutate(label_df = fct_inorder(label_df),
         label_gamma = fct_inorder(label_gamma),
         label_eta = fct_inorder(label_eta),
         label_method = fct_inorder(label_method))


#######################################
# combine results into tidy format.
#######################################

comb_res_slurm <- function(i){
  
  est <- NULL
  
  for(j in 1:nsim){
    
    if(file.exists(paste0("../scen", i, "/est", j,".rds"))){
      temp <- readRDS(paste0("../scen", i, "/est", j,".rds"))
      temp <- temp %>%
        mutate(isim = j, mod_type = "sim") 
      est <- rbind(est, temp)
    }
    
    
  }
  
  saveRDS(est, paste0("../scen", i, "_est.rds"))
  
  true <- mod_est(mod = true_mod, estimand = NULL, t_vec, rmst_vec, nonprop = FALSE,
                  return_model = TRUE, save_file = NULL)
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("../scen", i, "_res.rds"))
  
}

nscen <- nrow(scenarios)
pars_comb <- tibble(i = 1:nscen)

sjob_comb_results <- slurm_apply(comb_res_slurm, 
                                 pars_comb, 
                                 jobname = paste0(jobname,"_comb"),
                                 nodes = 20, cpus_per_node = 4, 
                                 submit = TRUE,
                                 global_objects = c( "mod_est", "create_res",
                                                     "true_mod", "nsim", "t_vec", "rmst_vec"),
                                 pkgs = c("dplyr", "tidyr", "ggplot2", "readr",
                                          "flexsurv", "survextrap"),
                                 slurm_options = list(time='02:00:00',
                                                      partition='core',
                                                      "mem-per-cpu"= '24G'))

system(check_status)

test <- readRDS("scen1_res.rds")
test <- readRDS( paste0("_rslurm_",jobname,"_comb/results_0.RDS"))
test
#View(test)
#summary(as.factor(test$estimand))
