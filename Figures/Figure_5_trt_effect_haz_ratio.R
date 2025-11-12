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

store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"

source("Functions/estimands.R")
source("Functions/performance.R")

jobname <- "cetux_trt3"
store_res <- paste0(store_directory, "simsurvextrap_trt3/") 
setwd(store_res)

####################################################
# Plot hazard ratios for non-PH models.
####################################################

#setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))

scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
nscen <- nrow(scenarios)

scen_df <- scenarios %>%
  filter(mspline_df %in% c(6,10),
         prior_hrsd_rate %in% c(1,5),
        mspline_bsmooth %in% TRUE,
         nonprop == TRUE,
         stan_fit_method == "mcmc")


store_plot_res <- paste0("plots/","plot_res_", jobname, ".rds")
res <- readRDS(store_plot_res)

res <- res %>%
  filter(isim <= 50) %>%
  filter(estimand == "hr") %>%
  filter(scenario_id %in% scen_df$scenario_id) %>%
  filter(t != 0)

length(unique(res$estimand_id))
length(unique(res$scenario_id))


for(smooth_model_type in c("random_walk", "exchangeable")) {
  #smooth_model_type <- "random_walk"
  if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
  if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
  
  
  scen_df <- scenarios %>%
    filter(smooth_model %in% smooth_model_type)
  
  #names(scen_df)
  #names(res)
  #head(res$label_hrsd)
  
  label_vec <- c("(a) Scenario 1: Constant effect (proportional hazards)",
                 "(b) Scenario 2: Waning effect", 
                 "(c) Scenario 3: Delayed then waning effect",
                 "(d) Scenario 4: Crossing survival curves")
  
  plot_hr <- res %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    mutate(label_scenario = case_when(
      hr_scenario == 1 ~ "Scenario 1: Constant effect \n (proportional hazards)",
      hr_scenario == 2 ~ "Scenario 2: Waning effect",
      hr_scenario == 3 ~ "Scenario 3: Delayed then waning effect",
      hr_scenario == 4 ~ "Scenario 4: Crossing survival curves"
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
    geom_segment(aes(x = 0, xend = 5.5, y = 1, yend = 1), colour = "gray55", linewidth = 0.65, linetype = "dashed")+
    geom_line()+
    scale_linewidth(range = c(0.4,1))+
    scale_alpha(range = c(0.1,1))+
    scale_color_manual(values = c("#830051", "black"))+
    facet_wrap(~label_scenario+label_method, ncol = 3, labeller = labeller(label_method = label_parsed))+
    scale_y_log10("Hazard ratio (log scale)", lim = c(0.42,2.25), 
                  breaks = c(0.50, 1.00, 2.00))+
    scale_x_continuous("Time (years)", breaks = 1:5)+
    guides(linewidth = "none", alpha = "none", colour = "none") 
  
  #plot_hr
  
  setwd(store_directory)
  
  pdf(file = paste0(store_directory, "figures/figure_trt_haz_",smooth_model_ext,".pdf"),   
      width = 5.2, 
      height = 6.4) 
  print(plot_hr)
  dev.off()
  
  
  jpeg(file = paste0(store_directory, "figures/figure_trt_haz_",smooth_model_ext,".jpg"),   
       width = 5.2, 
       height = 6.4,
       units = 'in',  
       res = 600)
  print(plot_hr)
  dev.off()
  
}

######################################

