#############################################################
# Load packages.
#############################################################

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

jobname <- "cetux_OS" # or "niv_PFS" 
store_res <- "directory/to/store/simulations"

#############################################################
# import trial data and call scripts to create functions.
#############################################################

source("Functions/simulate_dgm.R")
source("Functions/estimands.R")

cetux <- readRDS("Data/cetuximab_OS.rds") 
niv <- readRDS("Data/nivolumab_PFS.rds") 

dir.create(store_res)
setwd(store_res)

#########################################################
# DGM: import trial data
#########################################################

if(jobname == "cetux_OS"){
  
  # Cetuximab OS
  surv_df <- cetux
  
  # end of trial follow-up time for rmst computation
  maxT_data <- max(surv_df$years)
  
  k_true <-  3
  true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
  
} else if(jobname == "niv_PFS"){
  
  # Nivolumab PFS
  surv_df <- niv
  
  # end of trial follow-up time for rmst computation
  maxT_data <- max(surv_df$years)
  
  k_true <-  6
  true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
  
}


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

# Specify the priors and M-splines to vary in simulation

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

store_scenarios <- paste0(store_res,"scenarios_", jobname, ".csv") 
write_csv(scenarios, file = store_scenarios)
scenarios <- read_csv(store_scenarios)

# Number of simulation scenarios
nscen <- nrow(scenarios)

# Generate prior and mspline objects for survextrap in 
# each scenario

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
# Simulate dataset from dgm
########################################################


ndgm <- nrow(scenarios_dgm)

for(i in 1:ndgm){
  
  true_mod <- get(scenarios_dgm$dgm_mod[i])
  
  # Generate the simulated clinical trial data using sim_dgm function.
  # may take a few minutes if nsim large.
  
  assign(paste0("sim_data",i), sim_dgm(model_class = "flexsurv",
                                       true_mod = true_mod,
                                       N = scenarios_dgm$N[i],
                                       nsim = nsim, 
                                       maxT = scenarios_dgm$maxT[i],
                                       seed = scenarios_dgm$dgm_seed[i]) )
  
}


############################################
# Create directories to store results
############################################

# Folders to store results for each scenario.
for(i in 1:nscen){
  if(!dir.exists(paste0("scen", i))){
    dir.create(paste0("scen", i))
  }
}

#######################################################
# Create tibble/data.frame of parameters to call from
#######################################################

# Time points for hazard and survival estimates
t_vec <- seq(from = 0, to = maxT_data, length.out = 100)

# Time points for rmst estimates
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

store_scenarios_nsim <- paste0(store_res,"scenarios_nsim_", jobname, ".csv") 
write_csv(scenarios_nsim, file = store_scenarios_nsim)
# scenarios_nsim <- read_csv(store_scenarios_nsim)

pars <- scenarios_nsim %>%
  select(-c(scenario_id, dgm_id, num)) %>%
  relocate(sim_data_chr) %>%
  mutate(num = row_number()) %>%
  slice(sample(1:n())) 

store_parameters <- paste0(store_res,"parameters_", jobname, ".csv")
write_csv(pars, file = store_parameters)
# pars <- read_csv(store_parameters)

######################################
## HPC Slurm job for model fitting
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
                "catchToList", "throw_sampler_warnings",
                "diags", 
                "rmst_vec", "t_vec",
                "mod_est", "true_mod")

attach_pkgs <- c("dplyr", "tidyr", "ggplot2", "readr",
                 "flexsurv", "survextrap", "rstan")

slurm_opts <- list(time='23:00:00',
                   partition='core',
                   "mem-per-cpu"= '16G')
 
# Use slurm apply function.
sjob_fit_est <- slurm_apply(fit_est_slurm, 
                                pars_slurm, 
                                jobname = jobname,
                                nodes = 20, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                global_objects = attach_obj,
                                pkgs = attach_pkgs ,
                                slurm_options = slurm_opts )

store_sjob <- paste0(store_res, "sjob_", jobname, ".rds")
saveRDS(sjob_fit_est, file = store_sjob)
# sjob_fit_est <- readRDS(file = store_sjob)

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
