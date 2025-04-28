### Figure #


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

#######################################################

jobname <- "trt6" 

user <- Sys.info()["user"]
store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")

source("R/simulate_dgm_treatment_effect.R")
source("R/simulate_dgm.R")
source("R/functions_treatment_effect.R")
source("R/fit_model.R")
source("R/estimands.R")
source("R/visualise.R")
source("R/performance.R")




user <- Sys.info()["user"]
store_res <- paste0(store_directory, "simsurvextrap_slurm_", jobname, "/") 
setwd(store_res)

#setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))

##########################################
# Plot RMST at 5-y performance measures
##########################################

#dir.create("plots")

#################################
# Get truth.
#################################


scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
nscen <- nrow(scenarios)
nscen
scenarios

for(i in scenarios$scenario_id){

  temp <- readRDS(paste0("scen", i, "_res.rds"))
  
  temp <- temp %>% 
    filter(t == 5) %>%
    filter(estimand == "irmst") %>%
    mutate(estimand_id = "irmst1") %>%
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

saveRDS(res, "plots/res.RDS")

# true values.
res %>%
  filter(isim == 0) %>%
  filter(!duplicated(value))

##########################################################
# Performance measures.
##########################################################


est_id_choose <- "irmst1"
perform_res <- NULL
for(i in 1:nscen){
  temp_res <- sim_perform(res %>%
                            filter(scenario_id == i) %>%
                            filter(estimand_id == est_id_choose), 
                          estimands = c("irmst")) %>%
    mutate(scenario_id = i)
  
  if(i == 1){ 
    perform_res <- temp_res
  } else {
    perform_res <- rbind(perform_res, temp_res)
  }
  
}
saveRDS(perform_res, "plots/perform_res.RDS")
#head(perform_res)

all_res <- 
  perform_res %>%
  filter(stat %in% c( "nsim", "bias", "rbias","empse", "mse", "modelse", "cover")) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round))

saveRDS(all_res, "plots/all_res.RDS")

####################################
# Plot.
####################################

res <- readRDS("plots/res.RDS")
perform_res <- readRDS("plots/perform_res.RDS")
all_res <- readRDS("plots/all_res.RDS")

names(scenarios)
summary(as.factor(scenarios$prior_hsd_rate))
summary(as.factor(scenarios$prior_hrsd_rate))

for(smooth_model_type in c("random_walk", "exchangeable")){
  #smooth_model_type <- "exchangeable"
  
  if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
  if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
  
  
  
  scen_df <- scenarios %>%
    filter(mspline_df %in% c(3,6, 10),
           prior_hsd_rate %in% c(1,5,20),
           stan_fit_method %in% "mcmc",
           hscale_default %in% TRUE,
           smooth_model %in% smooth_model_type,
           mspline_bsmooth %in% TRUE) %>%
    arrange(mspline_df) %>%
    arrange(prior_hsd_rate) %>%
    arrange(hr_scenario) %>%
    mutate(label_hrsd = paste0("`G`*`amma`*`(`*", 
                               prior_hrsd_shape, "*`,`*" , 
                               prior_hrsd_rate, "*`)`")) %>%
    mutate(label_hrsd = case_when(
      is.na(prior_hrsd_rate) ~ "PH~model",
      .default = label_hrsd
    )) %>%
    select(scenario_id,
           hr_scenario,
           mspline_df,
           label_hrsd) %>%
    mutate(y_height = 0) %>%
    mutate(mspline_df = as.character(mspline_df)) %>%
    mutate(fontface_chr = "plain") %>%
    add_row(tibble_row(scenario_id = 0, 
                       mspline_df = "M-spline df    ", 
                       label_hrsd = paste0("bold(tau", "~", "prior)"), 
                       y_height = 1,
                       fontface_chr = "bold"))# %>%
  #  mutate(y_height = y_height + 7) 
  
  #names(scen_df)
  #View(scen_df)
  
  #########################################
  # Plot 1. Model
  #########################################
  names(scen_df)
  
  label_vec <- c("(a) Scenario 1: Constant effect (proportional hazards)",
                 "(b) Scenario 2: Waning effect", 
                 "(c) Scenario 3: Delayed then waning effect")
  
  for(scenario_num in 1:3){
    #scenario_num <- 1
    
    plot_df <- scen_df %>%
      filter(hr_scenario == scenario_num | is.na(hr_scenario)) %>%
      arrange(-y_height) %>%
      mutate(y_height = -row_number()) 
    
    y_limits <- c(min(plot_df$y_height), max(plot_df$y_height)+0.25)
    
    plot1 <- plot_df %>%
      ggplot(aes(x = 0, y = y_height)) +
      theme_classic()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.y=element_blank(),
            plot.margin = unit(c(0.1, 0, 0, 0.2), "cm"))+
      geom_text(aes(x= 0, y = y_height, label =  mspline_df, fontface = fontface_chr), size = 2.5)+
      geom_text(aes(x= 0.075, y = y_height, label =  label_hrsd, fontface = fontface_chr), parse = T,size = 2.5)+
      scale_x_continuous("", limits = c(-0.05, .15), expand = c(0, 0))+
      ylim(y_limits)
    plot1
    
    #########################################
    # Plot 2. Forest/lollipop with RMST
    #########################################
    
    est_mean <- res %>%
      filter(estimand_id == est_id_choose) %>%
      filter(isim != 0) %>%
      group_by(scenario_id) %>%
      summarise(estimand_mean = mean(value))
    
    head(est_mean)
    
    est_true <- res %>%
      filter(estimand_id == est_id_choose) %>%
      filter(isim == 0) %>%
      group_by(scenario_id) %>%
      summarise(estimand_true = mean(value)) %>%
      ungroup()
    
    
    head(est_true)
    
    bias <- all_res %>%
      filter(stat == "bias") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id)
    head(bias)
    
    plot_df2 <- cbind(est_mean, 
                      est_true %>% select(estimand_true),
                      bias %>% select(-c(scenario_id))) %>%
      right_join(plot_df, by = "scenario_id") %>%
      filter(!is.na(hr_scenario)) #%>%
    
    #plot_df$scenario_id
    
    #xlim1 <-
    #0.55-0.36
    xlim1 <- c(0.31, 0.50)
    xbreaks1 <- seq(from = 0.35, to = 0.45, by = 0.05)
    xlim2 <- c(0.36, 0.55)
    xbreaks2 <- seq(from = 0.35, to = 0.50, by = 0.05)
    xlim3 <- c(0.41, 0.60)
    xbreaks3 <- seq(from = 0.40, to = 0.55, by = 0.05)
    
    
    plot2 <- plot_df2 %>%
      ggplot(aes(x = estimand_mean, y = y_height)) +
      theme_classic()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x = element_text( size = 6, angle = 45, vjust = 0.5, hjust=0.5),
            axis.title.x = element_text( face="bold",size = 8),
            plot.margin = unit(c(0.1, 0, 0, 0), "cm")) + 
      geom_segment(x = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                   xend = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                   y = max(plot_df2$y_height)+0.5, 
                   yend = min(plot_df2$y_height)-0.5,
                   colour= "#830051")+
      geom_segment(aes(x= estimand_true, xend = estimand_mean, y = y_height, yend = y_height),
                   alpha = 0.4, colour= "#830051")+
      #  geom_errorbar(aes(xmin=estimand_true+est_low, xmax=estimand_true+est_high),
      #                width=1,                    # Width of the error bars
      #                position=position_dodge(.9),
      #                colour = "red",
      #                alpha = 0.7)+
      geom_point(aes(x=estimand_true+est_low, y = y_height), shape = 91, size = 2.5, colour= "#830051")+
      geom_point(aes(x=estimand_true+est_high, y = y_height), shape = 93, size = 2.5, colour= "#830051")+
      geom_point(alpha = 0.7, colour= "#830051") +
      #scale_x_continuous("Difference in \n RMST at 5-y", 
      #                   limits = c(0.408, 0.448),
      #                   breaks = seq(from = 0.41, to = 0.45, by = 0.01),
      #                   expand = c(0, 0))+
      scale_x_continuous("Difference in \n RMST at 5-y", 
                         limits = get(paste0("xlim", scenario_num)),
                         breaks = get(paste0("xbreaks", scenario_num)),
                         expand = c(0, 0))+
      ylim(y_limits)
    
    plot2
    
    #########################################
    # Plot 3. Forest/lollipop with Coverage
    #########################################
    
    bias_round <- all_res %>%
      filter(stat == "bias") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id) %>%
      mutate(est_round = sprintf("%.3f", round(est, 3)),
             mcse_round = sprintf("%.3f", round(mcse, 3)),
             est_round = as.character(est_round),
             mcse_round = as.character(mcse_round)) %>%                                                    
      mutate(bias_chr = paste0(est_round, " (", mcse_round, ")")) 
    
    relbias_round <- all_res %>%
      filter(stat == "rbias") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id) %>%
      mutate(est_round = sprintf("%.2f", round(100*est, 2)),
             mcse_round = sprintf("%.2f", round(100*mcse, 2)),
             est_round = as.character(est_round),
             mcse_round = as.character(mcse_round)) %>%
      mutate(relbias_chr = paste0(est_round, "%", " (", mcse_round, "%)"))
    
    
    
    empse_round <- perform_res %>%
      filter(stat == "empse") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id) %>%
      mutate(est_round = sprintf("%.3f", round(est, 3)),
             mcse_round = sprintf("%.3f", round(mcse, 3)),
             est_round = as.character(est_round),
             mcse_round = as.character(mcse_round)) %>%
      mutate(empse_chr = paste0(est_round, " (", mcse_round, ")"))
    
    modse_round <- perform_res %>%
      filter(stat == "modelse") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id) %>%
      mutate(est_round = sprintf("%.3f", round(est, 3)),
             mcse_round = sprintf("%.4f", round(mcse, 4)),
             est_round = as.character(est_round),
             mcse_round = as.character(mcse_round)) %>%
      mutate(modse_chr = paste0(est_round, " (", mcse_round, ")"))
    
    #head(modse_round)
    
    plot_df3 <- bias_round %>%
      select(scenario_id, bias_chr) %>%
      left_join(relbias_round %>% select(scenario_id, relbias_chr), join_by(scenario_id)) %>%
      left_join(empse_round %>% select(scenario_id, empse_chr), join_by(scenario_id)) %>%
      left_join(modse_round %>% select(scenario_id, modse_chr), join_by(scenario_id)) %>%
      right_join(plot_df, by = "scenario_id") %>%
      filter(!is.na(hr_scenario))
    
    # plot_df3 <- cbind( bias_round %>% select(scenario_id, bias_chr),
    #                    relbias_round %>% select(relbias_chr)) %>%
    #   right_join(plot_df, by = "scenario_id") %>%
    #   filter(!is.na(hr_scenario))
    
    plot3 <- plot_df3 %>%
      mutate(fontface_chr = "plain") %>%
      add_row(bias_chr = "Bias \n (MCSE)",
              relbias_chr = "Relative bias (%) \n (MCSE)",
              empse_chr = "Empirical SE \n (MCSE)",
              modse_chr = "Model SD \n (MCSE)",
              y_height = -1,
              fontface_chr = "bold") %>%
      ggplot(aes(x = 0, y = scenario_id)) +
      theme_classic()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.y=element_blank(),
            plot.margin = unit(c(0.1, 0, 0, 0), "cm") )+
      #  geom_text(aes(x= 0, y = y_height, label =  bias_chr, fontface = fontface_chr), size = 2.5)+
      #  geom_text(aes(x= 0.075, y = y_height, label =  relbias_chr, fontface = fontface_chr), size = 2.5)+
      #  scale_x_continuous("", limits=  c(-0.05, .12), expand = c(0, 0))+
      geom_text(aes(x= 0, y = y_height, label =  bias_chr, fontface = fontface_chr), size = 2.5)+
      geom_text(aes(x= 0.08, y = y_height, label =  relbias_chr, fontface = fontface_chr), size = 2.5)+
      geom_text(aes(x= 0.165, y = y_height, label =  empse_chr, fontface = fontface_chr), size = 2.5)+
      geom_text(aes(x= 0.24, y = y_height, label =  modse_chr, fontface = fontface_chr), size = 2.5)+
      scale_x_continuous("", limits=  c(-0.05, .285), expand = c(0, 0))+
      ylim(y_limits)
    plot3
    
    #########################################
    # Plot 4. Forest/lollipop with Coverage
    #########################################
    
    
    plot_df4 <- all_res %>%
      filter(stat == "cover") %>%
      mutate(est_high = est + 1*mcse,
             est_low = est - 1*mcse) %>%
      select(est, mcse, est_low, est_high, scenario_id) %>%
      mutate(est_true = 0.95) %>%
      right_join(plot_df, by = "scenario_id") %>%
      filter(!is.na(hr_scenario))
    
    #plot_df4
    plot4 <- plot_df4 %>%
      ggplot(aes(x = est, y = y_height)) +
      theme_classic()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5, hjust=0.5),
            axis.title.x = element_text(face="bold", size = 8),
            plot.margin = unit(c(0.1, 0, 0, 0), "cm")) +
      geom_segment(x = 0.95, 
                   xend = 0.95, 
                   y = max(plot_df2$y_height)+0.5, 
                   yend = min(plot_df2$y_height)-0.5, colour= "#830051")+
      geom_segment(aes(x= est_true, xend = est, y = y_height, yend = y_height ),
                   alpha = 0.4, colour= "#830051")+
      # geom_errorbar(aes(xmin=est_low, xmax=est_high),
      #                width=1,                    # Width of the error bars
      #                position=position_dodge(.9),
      #                colour = "black",
      #                alpha = 0.6)+
      geom_point(aes(x=est_low, y = y_height), shape = 91, size = 2.5, colour= "#830051")+
      geom_point(aes(x=est_high, y = y_height), shape = 93, size = 2.5, colour= "#830051")+
      geom_point(alpha = 0.7, colour= "#830051") +
      scale_x_continuous("Coverage", limits = c(0.80, 1.00), labels = scales::percent,
                         breaks = seq(from = 0.80, to = 1.00, by = 0.05),
                         expand = c(0, 0))+
      ylim(y_limits)
    plot4
    
    
    ##################################################
    ##################################################
    
    wid1 <- 2.6
    wid1null <- -0.4
    wid2 <- 2.0
    wid2null <- -0.2
    wid3 <- 5.4
    wid3null <- -0.2
    wid4 <- 2.0
    wid4null <- 0.2
    
    figure_title <- ggdraw() + 
      draw_label(
        label_vec[scenario_num],
        size=10, 
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        plot.margin = margin(0, 0, 0, 7)
      )
    
    
    figure_row <- plot_grid(plot1, NULL,plot2, NULL, plot3, NULL, plot4, NULL, nrow = 1,
                            rel_widths = c(wid1, wid1null, wid2, wid2null, wid3,
                                           wid3null, wid4, wid4null),
                            align = "h")
    
    figure <- plot_grid(figure_title, figure_row,
                        ncol = 1,
                        rel_heights = c(0.1, 1))
    
    assign(paste0("figure", scenario_num), figure)
    saveRDS(figure, paste0("plots/figure2_",scenario_num,"_", smooth_model_ext ,".RDS"))
  }
  
  figure1
  figure2
  figure3
  
  #figure_cetux <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", "master_cetux10", "/plots/figure1.RDS"))
  #figure_niv <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", "master_niv3", "/plots/figure1.RDS"))
  
  figure_all <- plot_grid(NULL, figure1, NULL, figure2, NULL, figure3,
                          rel_heights = c(0.01, 0.4,0.01,  0.4, 0.01, 0.4), 
                          ncol = 1)
  
  pdf(file = paste0(store_directory, "figures/figure_trt_perform_",smooth_model_ext,".pdf"),   
      width = 7.5, 
      height = 8) 
  print(figure_all)
  dev.off()
  
  tiff(file = paste0(store_directory, "figures/figure_trt_perform_",smooth_model_ext,".tiff"),   
       width = 7.5, 
       height = 8,
       units = 'in',  
       res = 1200, 
       compression = "lzw")
  print(figure_all)
  dev.off()
  
}
