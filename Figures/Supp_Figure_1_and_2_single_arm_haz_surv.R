#Supp

###################

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


cetux_jobname <- "cetux16"
niv_jobname <- "niv8" 

scen_name <- "cetux"

jobname <- get(paste0(scen_name, "_jobname"))

user <- Sys.info()["user"]
store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")

source("R/fit_model.R")
source("R/estimands.R")
source("R/visualise.R")
source("R/performance.R")

##########################################
# Plot hazard and survival curves.
##########################################

setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))

#dir.create("plots")

#################################
# Get truth.
#################################

if(scen_name == "cetux"){
  # From Cetux control arm
  surv_df <- cetux[cetux$treat=="Cetuximab",]
  
  # end of trial follow-up time for rmst computation
  
  maxT_data <- max(surv_df$years)
  k_true <-  3
  true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
  t_vec <- seq(from = 0, to = maxT_data, length.out = 100)
  rmst_vec <- sort(c(maxT_data, seq(from = 2, to = 6, by = 1)))
}

if(scen_name == "niv"){
  surv_df <- read_csv(here::here("nivo_data/nivolumab_digiKM_pfs.csv"),
                      show_col_types = FALSE)  %>%
    filter(treat==2) %>%
    mutate(years = time/12) %>%
    rename(d = status)
  
  # end of trial follow-up time for rmst computation
  
  maxT_data <- max(surv_df$years)
  k_true <-  6
  true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
  t_vec <-sort(c(seq(from = 0, to = maxT_data, length.out = 100),
                 seq(from = 0.14, to = 0.26, length.out = 1e4 )))
  rmst_vec <- sort(c(maxT_data, seq(from = 2, to = 6, by = 1)))
}


true <- mod_est(mod = true_mod, estimand = NULL, t_vec, rmst_vec, nonprop = FALSE,
                return_model = TRUE, save_file = NULL)

true <- true %>%
  mutate(isim = 0, mod_type = "true")

#########################################
# Combine results across simulations.
#########################################

scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
nscen <- nrow(scenarios)

nsim <- 50

for(i in 1:nscen){
  
  est <- NULL
  
  for(j in 1:nsim){
    # print progress track every 10th iteration
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
    if(file.exists(paste0("scen", i, "/est", j,".rds"))){
      #i <- j <- 1
      temp <- readRDS(paste0("scen", i, "/est", j,".rds"))
      
      temp <- temp %>%
        mutate(isim = j, mod_type = "sim") 
      
      est <-  rbindlist(list(est, temp))
    }
  }
  est <- data.frame(est)
  saveRDS(est, paste0("scen", i, "_est_sample.rds"))
  
  # append true values alongside estimates.
  res <- rbind(true, est) %>%
    group_by(estimand, isim) %>%
    mutate(id = row_number(),
           estimand_id = paste0(estimand, id)) %>%
    relocate(estimand_id) %>%
    select(-id) %>%
    ungroup()
  
  saveRDS(res, paste0("scen", i, "_res_sample.rds"))
}

#}
#View(res)

res_test <- readRDS("scen6_est_sample.rds")
res_test
#View(res_test)
#summary(as.factor(res_test$estimand_id))

###################################################
# Combine data together.
###################################################

for(i in scenarios$scenario_id){
  
  temp <- readRDS(paste0("scen", i, "_res_sample.rds"))
  
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
res



####################################
# Plot.
####################################

names(scenarios)
summary(as.factor(scenarios$prior_hsd_rate))
names(scenarios)


smooth_model_type <- "random_walk"
#smooth_model_type <- "exchangeable"

if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"


scen_df <- scenarios %>%
  filter(mspline_df %in% c(3,6, 10),
         prior_hsd_rate %in% c(1,5,20),
         stan_fit_method %in% "mcmc",
         smooth_model %in% smooth_model_type,
         hscale_default %in% TRUE,
         mspline_bsmooth %in% TRUE) %>%
  mutate(label_df = paste0("`df=`*", mspline_df)) %>%
  mutate(label_hsd = paste0("sigma*","`~`*`G`*`amma`*`(`*", 
                            prior_hsd_shape, "*`,`*" , 
                            prior_hsd_rate, "*`)`")) %>% 
  mutate(label_scen = paste0("`df=`*", mspline_df, "*`,`*", "~sigma*","`~`*`G`*`amma`*`(`*", 
                             prior_hsd_shape, "*`,`*" , 
                             prior_hsd_rate, "*`)`" )) %>%
  select(scenario_id,
         mspline_df,
         prior_hsd_rate,
         label_df,
        label_hsd,
        label_scen,
        prior_hsd_rate)

#names(scen_df)
#View(scen_df)

#########################################
# Plot 1. Model
#########################################

#head(res)
head(res)  
head(scen_df)

summary(as.factor(res$scenario_id))
  
survplot <- res %>%
    filter(estimand == "survival") %>%
    filter(t != 0) %>%
    filter(isim <= 50) %>% # sample of 50 iterations.
    left_join(scen_df, by = "scenario_id") %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    arrange(prior_hsd_rate) %>%
    arrange(mspline_df) %>%
    mutate(label_df = fct_inorder(label_df)) %>%
    mutate(label_hsd = fct_inorder(label_hsd)) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
  	mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
  	mutate(line_colour = if_else(isim  == 0, 1, 0),
         line_colour = as.factor(line_colour)) %>%
    ggplot(aes(x = t, y= value, alpha = line_alpha, colour = line_colour, 
		group = isim, linewidth = line_width))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
    geom_line()+
    scale_linewidth(range = c(0.4,1))+
    scale_alpha(range = c(0.1,1))+
  scale_color_manual(values = c("#830051", "black"))+
    facet_wrap(~label_df+label_hsd, ncol = 3, labeller = label_parsed)+
    scale_x_continuous("Time (years)", breaks = 1:5)+
    scale_y_continuous("Survival", limit = c(0,1), labels = scales::percent)+
    guides(linewidth = "none", alpha = "none", colour = "none")  
  
if(scen_name == "cetux"){ 
  y_limits <-  c(0, 0.7)
  y_breaks <- seq(0, 0.75, by = 0.2)
}
if(scen_name == "niv"){
  y_limits <-  c(0, 3)
  y_breaks <- seq(0, 3, by = 1)
}

hazplot <- res %>%
  filter(estimand == "hazard") %>%
  filter(t != 0) %>%
  filter(isim <= 50) %>% # sample of 50 iterations.
  left_join(scen_df, by = "scenario_id") %>%
  filter(scenario_id %in% scen_df$scenario_id) %>%
  arrange(prior_hsd_rate) %>%
  arrange(mspline_df) %>%
  mutate(label_df = fct_inorder(label_df)) %>%
  mutate(label_hsd = fct_inorder(label_hsd)) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
  	mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
  	mutate(line_colour = if_else(isim  == 0, 1, 0),
         line_colour = as.factor(line_colour)) %>%
  ggplot(aes(x = t, y= value, alpha = line_alpha, colour = line_colour, 
		group = isim, linewidth = line_width))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 6)) +
  geom_line()+
  scale_linewidth(range = c(0.4,1))+
  scale_alpha(range = c(0.1,1))+
  scale_color_manual(values = c("#830051", "black"))+
  facet_wrap(~label_df+label_hsd, ncol = 3, labeller = label_parsed)+
  scale_x_continuous("Time (years)", breaks = 1:5)+
  scale_y_continuous("Hazard",limits = y_limits, breaks = y_breaks)+
  guides(linewidth = "none", alpha = "none", colour = "none")  

figure_all <- plot_grid(NULL, survplot, NULL, hazplot, 
                        rel_heights = c(0.02, 0.5, 0.02, 0.5),
                        labels = c("(a) Survival", "", "(b) Hazard", ""),
                        label_size = 10,
                        ncol = 1,
                        align = "v")

figure_all
#getwd()


saveRDS(figure_all, paste0("plots/figure_haz_surv_",smooth_model_ext,".RDS"))

#####################################


pdf(file = paste0(store_directory, "figures/figure_haz_surv_",scen_name,"_",smooth_model_ext,".pdf"),   
    width = 6.5, 
    height = 9) 

figure_all

dev.off()

tiff(file = paste0(store_directory, "figures/figure_haz_surv_",scen_name,"_",smooth_model_ext,".tiff"),   
    width = 6.5, 
    height = 9,
     units = 'in',  
     res = 1200, 
     compression = "lzw")

figure_all

dev.off()


