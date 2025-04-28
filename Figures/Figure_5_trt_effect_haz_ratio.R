#############################################
# Hazard ratio plots, across tau priors.
#############################################

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


#######################################################

jobname <- "trt6" 

user <- Sys.info()["user"]
store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")

####################################################
# Plot hazard ratios for non-PH models.
####################################################

setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))

scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
nscen <- nrow(scenarios)

store_plot_res <- paste0("plots/","plot_res_", jobname, ".rds")
res <- readRDS(store_plot_res)


for(smooth_model_type in c("random_walk", "exchangeable")) {
  
  if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
  if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
  
  
  scen_df <- scenarios %>%
    filter(mspline_df %in% c(6,10),
           prior_hrsd_rate %in% c(1,5),
           smooth_model %in% smooth_model_type,
           mspline_bsmooth %in% TRUE)
  
  names(scen_df)
  names(res)
  head(res$label_hrsd)
  
  label_vec <- c("(a) Scenario 1: Constant effect (proportional hazards)",
                 "(b) Scenario 2: Waning effect", 
                 "(c) Scenario 3: Delayed then waning effect")
  
  plot_hr <- res %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    filter(estimand == "hr") %>%
    filter(nonprop == TRUE) %>%
    filter(t != 0) %>%
    filter(isim <= 50) %>%
    mutate(label_scenario = case_when(
      hr_scenario == 1 ~ "Scenario 1: Constant effect \n (proportional hazards)",
      hr_scenario == 2 ~ "Scenario 2: Waning effect",
      hr_scenario == 3 ~ "Scenario 3: Delayed then waning effect"
    )) %>% 
    #  mutate(hr_scenario = paste0("Scenario~", hr_scenario)) %>%
    mutate(label_method = paste0("df*`=`*", mspline_df ,"*`,`*", "~tau*`~`*" ,label_hrsd)) %>%
    arrange(hr_scenario, prior_hrsd_rate) %>%
    # mutate(label_all = fct_inorder(label_all)) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_colour = if_else(isim  == 0, 1, 0),
           line_colour = as.factor(line_colour)) %>%
    mutate(isim = as.factor(isim))  %>%
    ggplot(aes(x = t, y= value, alpha = line_alpha, colour = line_colour , group = isim, linewidth = line_width))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 5, 
                                      margin = margin( b = 3, t = 3))) +
    geom_line()+
    scale_linewidth(range = c(0.4,1))+
    scale_alpha(range = c(0.1,1))+
    scale_color_manual(values = c("#830051", "black"))+
    facet_wrap(~label_scenario+label_method, ncol = 3, labeller = labeller(label_method = label_parsed))+
    scale_y_continuous("Hazard ratio",limits = c(0.2, 1.4), breaks = seq(0.25, 1.5, by = 0.25))+
    scale_x_continuous("Time (years)", breaks = 1:5)+
    guides(linewidth = "none", alpha = "none", colour = "none") 
  
  plot_hr
  
  setwd(store_directory)
  
  pdf(file = paste0(store_directory, "figures/figure_trt_haz_",smooth_model_ext,".pdf"),   
      width = 5, 
      height = 5) 
  print(plot_hr)
  dev.off()
  
  
  tiff(file = paste0(store_directory, "figures/figure_trt_haz_",smooth_model_ext,".tiff"),   
       width = 5, 
       height = 5,
       units = 'in',  
       res = 1200, 
       compression = "lzw")
  print(plot_hr)
  dev.off()
  
}

######################################

