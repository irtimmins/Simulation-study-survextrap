#####################################################
library(rsimsum)
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

jobname <- "trt5_rw" 

user <- Sys.info()["user"]
store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"

source("Functions/estimands.R")
source("Functions/performance.R")

user <- Sys.info()["user"]
store_res <- paste0(store_directory, "simsurvextrap_slurm_", jobname, "/") 
setwd(store_res)

#####################################################

setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))

scenarios_flex <- expand_grid(dgm = 1:3,
                              k_knots = c(0,1,2,4,8)) %>%
  mutate(flex_scenario_id = row_number()) %>%
  relocate(flex_scenario_id) %>%
  mutate(df = k_knots + 2,
         df_tvc = df,
         new_id = paste0(flex_scenario_id, "_", "flex")) %>%
  select(dgm, df, df_tvc, new_id)

#head(scenarios_flex)
scenarios_rstpm <- expand_grid(dgm = 1:3,
                               df_base = c(1,2,3,5,9),
                               df_tvc = c(1,2,3,5,9)) %>%
  mutate(rstpm2_scenario_id = row_number()) %>%
  relocate(rstpm2_scenario_id) %>%
  mutate(df = df_base + 1,
         df_tvc = df_tvc + 1,
         new_id = paste0(rstpm2_scenario_id, "_rstpm2")) %>%
  select(dgm, df, df_tvc, new_id)

head(scenarios_rstpm)

head(scenarios_flex)
scen_df <- rbind(scenarios_flex %>%  mutate(model = "flexsurv"), 
                 scenarios_rstpm %>% mutate(model = "rstpm2")) %>%
  filter(df >= df_tvc)

#View(scen_df)
#View(scen_df)

################################################################
# Read in res files.
################################################################


for(i in 1:nrow(scenarios_flex)){
 # i<- 1
  
  temp <- readRDS(paste0("scen_flex", i, "_res_rmst.rds"))
  temp <- temp %>%
    mutate(new_id = paste0(i, "_flex"))

  #head(temp)

  if(i == 1) res_flex <- temp
  if(i > 1) res_flex <- rbind(res_flex, temp)
  print(i)

}

#head(res_flex)

for(i in 1:nrow(scenarios_rstpm)){
  
  temp <- readRDS(paste0("scen_rstpm2_", i, "_res_rmst.rds"))
  temp <- temp %>%
    mutate(new_id = paste0(i, "_rstpm2"))
  
  head(temp)
  
  if(i == 1) res_rstpm2 <- temp
  if(i > 1) res_rstpm2 <- rbind(res_rstpm2, temp)
  print(i)
  
}

res <- rbind(res_flex, res_rstpm2) %>%
    filter(!is.na(value)) %>%
  filter(value != 0) %>%
  filter(!is.na(value_se) | isim == 0)

##########################################################
# Performance measures.
##########################################################

est_id_choose <- "irmst1"
perform_res <- NULL

for(i in 1:nrow(scen_df)){
  temp_res <- sim_perform(res %>%
                            filter(new_id == scen_df$new_id[i]) %>%
                            filter(estimand_id == est_id_choose), 
                          estimands = c("irmst")) %>%
    mutate(new_id = scen_df$new_id[i])
  
  if(i == 1){ 
    perform_res <- temp_res
  } else {
    perform_res <- rbind(perform_res, temp_res)
  }
  
}

#head(perform_res)
all_res <- 
  perform_res %>%
  filter(stat %in% c( "nsim", "bias", "rbias","empse", "mse", "modelse", "cover")) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round))

##################################################################

label_vec <- c("(a) Scenario 1: Constant effect (proportional hazards)",
               "(b) Scenario 2: Waning effect", 
               "(c) Scenario 3: Delayed then waning effect")



for(scenario_num in 1:3){
  #scenario_num <- 1
  #head(scen_df)
  
  plot_df <- scen_df %>%
    filter(dgm == scenario_num) %>%
    arrange(desc(model)) %>%
    arrange(df_tvc) %>%
    arrange(df) %>%
    mutate(df = as.character(df),
           df_tvc = as.character(df_tvc)) %>%
    select(new_id, model, df, df_tvc) %>%
    mutate(y_height = -row_number()) %>%
    mutate(fontface_chr = "plain") %>%
    add_row(tibble_row(new_id = "0", 
                       model = "R Package", 
                       df = "Baseline \n spline, df", 
                       df_tvc = "Time-varying \n effect spline, df", 
                       y_height = 1,
                       fontface_chr = "bold")) 
  
  y_max <- max(plot_df %>%
                 select(y_height) %>%
                 pull())+0.5
  
  y_min <- min(plot_df %>%
                 select(y_height) %>%
                 pull())
  
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
    geom_text(aes(x= 0.02, y = y_height, label =  model, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.16, y = y_height, label =  df, fontface = fontface_chr), parse = F,size = 2.5)+
    geom_text(aes(x= 0.34, y = y_height, label =  df_tvc, fontface = fontface_chr), parse = F,size = 2.5)+
        scale_x_continuous("", limits = c(-0.05, .495), expand = c(0, 0))+
    ylim(c(y_min,y_max))
  
  plot1
  
  #########################################
  # Plot 2. Forest/lollipop with RMST
  #########################################
  
  #head(res)
  #summary(as.factor(res$estimand_id))
  est_mean <- res %>%
    filter(isim != 0) %>%
    group_by(new_id) %>%
    summarise(estimand_mean = mean(value)) %>%
    ungroup()
  
  head(est_mean)
  
  est_true <- res %>%
    filter(isim == 0) %>%
    group_by(new_id) %>%
    summarise(estimand_true = mean(value)) %>%
    ungroup()
  
  
  head(est_true)
  
  bias <- all_res %>%
    filter(stat == "bias") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, new_id)
  head(bias)
  
  plot_df2 <-  est_mean %>%
    left_join(est_true, join_by(new_id)) %>%
    left_join(bias, join_by(new_id)) %>%
    right_join(plot_df, join_by(new_id)) 
  
  head(plot_df2)
  #plot_df$scenario_id
  
  #xlim1 <-
  xlim1 <- c(0.31, 0.55)
  xbreaks1 <- seq(from = 0.35, to = 0.50, by = 0.05)
  xlim2 <- c(0.36, 0.60)
  xbreaks2 <- seq(from = 0.40, to = 0.55, by = 0.05)
  xlim3 <- c(0.36, 0.60)
  xbreaks3 <- seq(from = 0.40, to = 0.55, by = 0.05)
  
  
  plot2 <- plot_df2 %>%
    filter(!is.na(estimand_mean)) %>%
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
                 y = max(plot_df2$y_height)-0.5, 
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
    ylim(c(y_min,y_max))
  
  
  plot2
  
  #########################################
  # Plot 3. Forest/lollipop with Coverage
  #########################################
  bias_round <- perform_res %>%
    filter(stat == "bias") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, new_id) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round)) %>%                                                    
    mutate(bias_chr = paste0(est_round, " (", mcse_round, ")")) 
  
  relbias_round <- perform_res %>%
    filter(stat == "rbias") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, new_id) %>%
    mutate(est_round = sprintf("%.2f", round(100*est, 2)),
           mcse_round = sprintf("%.2f", round(100*mcse, 2)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round)) %>%
    mutate(relbias_chr = paste0(est_round, "%", " (", mcse_round, "%)"))
  
  
empse_round <- perform_res %>%
  filter(stat == "empse") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%
  mutate(empse_chr = paste0(est_round, " (", mcse_round, ")"))

modse_round <- perform_res %>%
  filter(stat == "modelse") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.4f", round(mcse, 4)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%
  mutate(modse_chr = paste0(est_round, " (", mcse_round, ")"))

#head(modse_round)

plot_df3 <- bias_round %>%
  select(new_id, bias_chr) %>%
  left_join(relbias_round %>% select(new_id, relbias_chr), join_by(new_id)) %>%
  left_join(empse_round %>% select(new_id, empse_chr), join_by(new_id)) %>%
  left_join(modse_round %>% select(new_id, modse_chr), join_by(new_id)) 
  
  plot3 <- plot_df3 %>%
    left_join(plot_df %>% select(new_id, y_height), join_by(new_id) ) %>%
    filter(!is.na(y_height)) %>%
    mutate(fontface_chr = "plain") %>%
    add_row(bias_chr = "Bias \n (MCSE)",
            relbias_chr = "Relative bias (%) \n (MCSE)",
            empse_chr = "Empirical SE \n (MCSE)",
            modse_chr = "Model SE \n (MCSE)",
            y_height = 1,
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
    geom_text(aes(x= 0, y = y_height, label =  bias_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.075, y = y_height, label =  relbias_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.160, y = y_height, label =  empse_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.235, y = y_height, label =  modse_chr, fontface = fontface_chr), size = 2.5)+
    scale_x_continuous("", limits=  c(-0.05, .28), expand = c(0, 0))+
    ylim(c(y_min,y_max))
  
  
  plot3
  
  #########################################
  # Plot 4. Forest/lollipop with Coverage
  #########################################
  
  
  plot_df4 <- perform_res %>%
    filter(stat == "cover") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, new_id) %>%
    mutate(est_true = 0.95) %>%
    right_join(plot_df, by = "new_id") 
  
  #plot_df4
  plot4 <- plot_df4 %>%
    filter(!is.na(y_height),
           !is.na(est)) %>%
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
                 y = max(plot_df2$y_height)-0.5, 
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
    scale_x_continuous("Coverage", limits = c(0.75, 1.00), labels = scales::percent,
                       breaks = seq(from = 0.75, to = 1.00, by = 0.05),
                       expand = c(0, 0))+
    ylim(c(y_min,y_max))
  
  plot4
  
  ##################################################
  ##################################################
  
  wid1 <- 3.7
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
                      rel_heights = c(0.08, 1))
  
  
  assign(paste0("figure", scenario_num), figure)
  saveRDS(figure, paste0("plots/supfigure_trt_bench_",scenario_num,".RDS"))
}

figure1
figure2
figure3

figure_all <- plot_grid(NULL, figure1, NULL, figure2, NULL, figure3,
                        rel_heights = c(0.01, 0.4,0.01,  0.4, 0.01, 0.4),
                        ncol = 1)
figure_all

pdf(file = paste0(store_directory, "figures/figure_trt_bench_perform.pdf"),   
    width = 8.5, 
    height = 12) 
figure_all
dev.off()

tiff(file = paste0(store_directory, "figures/figure_trt_bench_perform.tiff"),   
     width = 8.5, 
     height = 12,
     units = 'in',  
     res = 1200, 
     compression = "lzw")

figure_all
dev.off()

