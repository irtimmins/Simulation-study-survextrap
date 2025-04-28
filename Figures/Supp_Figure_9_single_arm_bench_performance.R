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

source("R/fit_model.R")
source("R/estimands.R")
source("R/visualise.R")
source("R/performance.R")

##########################################
# Combine data for flexsurv.
##########################################

#for(jobname in c("master_cetux12", "master_niv4_smooth" )){
jobname <- "cetux16"
#jobname <- "niv8" 
  
  setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))
  
  
  nsim <- 1000
  
  scenarios_flex <- expand_grid(k_knots = c(0,1,2,4,8)) %>%
    mutate(flex_scenario_id = row_number()) %>%
    relocate(flex_scenario_id)
  
  nscen_flex <- nrow(scenarios_flex)

   for(i in 1:nscen_flex){
    est <- NULL
    
    for(j in 1:nsim){
      # print progress track every 10th iteration
      if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
      if(file.exists(paste0("scen_flex", i, "/flex_est", j,".rds"))){
        temp <- readRDS(paste0("scen_flex", i, "/flex_est", j,".rds"))
        temp <- temp %>%
          mutate(isim = j, mod_type = "sim") 
        est <-  rbindlist(list(est, temp))
      }
    }
    est <- data.frame(est) %>%
      select(-trt)
    
    saveRDS(est, paste0("flex_scen", i, "_est.rds"))
   
    # read in true values of estimands.
    true <- readRDS("true_est.RDS")
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
    
    print(nrow(est))
    print(nrow(res))
    saveRDS(res, paste0("flex_scen", i, "_res.rds"))
    
  }
  
  #}


##########################################
# Combine data for rstpm2.
##########################################

  nsim <- 1000
  
  scenarios_rstpm <- expand_grid(df = c(1,2,3,5,9)) %>%
    mutate(rstpm2_scenario_id = row_number()) %>%
    relocate(rstpm2_scenario_id)
  
  
  nscen_rstpm <- nrow(scenarios_rstpm)
  
  for(i in 1:nscen_rstpm){
    #i <- 5
    est <- NULL
    
    for(j in 1:nsim){
      # print progress track every 10th iteration
      if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
      if(file.exists(paste0("scen_rstpm2_", i, "/rstpm2_est", j,".rds"))){
        temp <- readRDS(paste0("scen_rstpm2_", i, "/rstpm2_est", j,".rds"))
        temp <- temp %>%
          mutate(isim = j, mod_type = "sim") 
        est <-  rbindlist(list(est, temp))
      }
    }
    est <- data.frame(est) %>%
      select(-trt) %>%
      mutate(value_se = (value_ci_high-value_ci_low)/(2*1.96))

    
    saveRDS(est, paste0("rstpm2_scen", i, "_est.rds"))
    
    # read in true values of estimands.
    true <- readRDS("true_est.RDS")
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
    
    saveRDS(res, paste0("rstpm2_scen", i, "_res.rds"))
    
  }
  

##########################################
# Plot RMST at 5-y performance measures
##########################################

cetux_jobname <-  "cetux16"
niv_jobname <- "niv8" 
  
for(jobname in c(cetux_jobname, niv_jobname )){
  

  setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))
  
  scenarios_flex
  scenarios_rstpm

  ###################################################
  # Combine data together.
  ###################################################
  
  for(i in scenarios_flex$flex_scenario_id){

    temp <- readRDS(paste0("flex_scen", i, "_res.rds"))
    #head(temp)
    temp <- temp %>% 
      mutate(flex_scenario_id = i)  %>%
      left_join(scenarios_flex, 
                by = join_by(flex_scenario_id)) %>%
      mutate(method = "flexsurv") %>%
      mutate(df = k_knots + 2) %>%
      select(-c(flex_scenario_id,k_knots))
   
     #View(test)
    
    if(i == scenarios_flex$flex_scenario_id[1]){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("scenario ", i, "/", max(scenarios_flex$flex_scenario_id)))
    res <- data.frame(res)
    res_flex <- res %>%
      filter(!is.na(value))
  }
  #head(res_flex)
  
  for(i in scenarios_rstpm$rstpm2_scenario_id){

    temp <- readRDS(paste0("rstpm2_scen", i, "_res.rds"))
        #head(temp)
    
    temp <- temp %>% 
      mutate(rstpm2_scenario_id = i)  %>%
      left_join(scenarios_rstpm, 
                by = join_by(rstpm2_scenario_id)) %>%
      mutate(method = "rstpm2") %>%
      mutate(df = df + 1) %>%
      select(-c(rstpm2_scenario_id))
    
    #View(test)
    
    if(i == scenarios_rstpm$rstpm2_scenario_id[1]){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("scenario ", i, "/", max(scenarios_rstpm$rstpm2_scenario_id)))
    res <- data.frame(res)
    res_rstpm2 <- res %>%
      filter(!is.na(value))
  }

  # head(res)
  
  res <- rbind(res_flex,
               res_rstpm2) %>% 
    filter(value < 4.5) %>% # deal with delta method error.
    group_by(df, method) %>%
    mutate(scenario_id = cur_group_id())
  
  
  #summary(as.factor(res$scenario_id))
  #summary(res)

  est_id_choose <- "rmst1"
  perform_res <- NULL
  for(i in 1:length(unique(res$scenario_id))){
    #i <- 1
    
    temp_res <- sim_perform(res %>%
                              filter(scenario_id == i) %>%
                              filter(estimand_id == est_id_choose), 
                            estimands = c("rmst")) %>%
      mutate(scenario_id = i)
    
    if(i == 1){ 
      perform_res <- temp_res
    } else {
      perform_res <- rbind(perform_res, temp_res)
    }
    
  }
  
  head(perform_res)
  #View(perform_res)
  all_res <- 
    perform_res %>%
    filter(stat %in% c( "nsim", "bias", "rbias","empse", "mse", "modelse", "cover")) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round))
  

  ####################################
  # Plot.
  ####################################
  

scen_df <- res %>%
  filter(!duplicated(scenario_id)) %>%
    select(scenario_id,
           df,
           method) %>%
    arrange(df) %>%
    mutate(df = as.character(df)) %>%
    ungroup() %>%
    mutate(y_height = -row_number()) %>%
    mutate(fontface_chr = "plain") %>%
    add_row(tibble_row(scenario_id = 0, 
                       df = "Spline df", 
                       method = "R package", 
                       y_height = 0.5,
                       fontface_chr = "bold")) 
  
y_max <- max(scen_df %>%
  select(y_height) %>%
  pull())

y_min <- min(scen_df %>%
               select(y_height) %>%
               pull())
#names(scen_df)
  #View(scen_df)
  
  #########################################
  # Plot 1. Model
  #########################################
  
  plot1 <- scen_df %>%
    ggplot(aes(x = 0, y = y_height)) +
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.y=element_blank())+
    geom_text(aes(x= 0, y = y_height, label =  method, fontface = fontface_chr), parse = F, size = 2.5)+
    geom_text(aes(x= 0.075, y = y_height, label =  df, fontface = fontface_chr), parse = F,size = 2.5)+
    scale_x_continuous("", limits = c(-0.05, .15), expand = c(0, 0))+
    ylim(c(y_min,y_max))
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
                    bias %>% select(-c(scenario_id)))
  
  head(plot_df2)
  # View(plot_df2)
  test <- plot_df2 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id))
  #View(test)
  
  if(jobname == cetux_jobname) {
    xlim <- c(3.10, 3.34) 
    xbreaks <- seq(from = 3.10, to = 3.34, by = 0.05)                       
  } else{
    xlim <- c(1.85, 2.14)
    xbreaks <- seq(from = 1.85, to = 2.14, by = 0.05)
  }
  plot2 <- plot_df2 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id)) %>%
    ggplot(aes(x = estimand_mean, y = y_height)) +
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text( size = 6, angle = 45, vjust = 0.5, hjust=0.5),
          axis.title.x = element_text( face="bold",size = 8)) + 
    geom_segment(x = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                 xend = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                 y = y_min-0.5, 
                 yend = y_max-1, colour = "#830051")+
    geom_segment(aes(x= estimand_true, xend = estimand_mean, y = y_height, yend = y_height),
                 alpha = 0.4, colour = "#830051")+
    #  geom_errorbar(aes(xmin=estimand_true+est_low, xmax=estimand_true+est_high),
    #                width=1,                    # Width of the error bars
    #                position=position_dodge(.9),
    #                colour = "red",
    #                alpha = 0.7)+
    geom_point(aes(x=estimand_true+est_low, y = y_height), shape = 91, size = 2.5, colour = "#830051")+
    geom_point(aes(x=estimand_true+est_high, y = y_height), shape = 93, size = 2.5, colour = "#830051")+
    geom_point(alpha = 0.7, colour= "#830051") +
    scale_x_continuous("RMST at 5-y", limits = xlim,
                       breaks = xbreaks,
                       expand = c(0, 0))+
    ylim(c(y_min,y_max))
  
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
  
  
  #summary(as.factor(all_res$stat))
  empse_round <- all_res %>%
    filter(stat == "empse") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, scenario_id) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round)) %>%
    mutate(empse_chr = paste0(est_round, " (", mcse_round, ")"))
  
  modse_round <- all_res %>%
    filter(stat == "modelse") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, scenario_id) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.4f", round(mcse, 4)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round)) %>%
    mutate(modse_chr = paste0(est_round, " (", mcse_round, ")"))
  
  
  
  plot_df3 <- cbind( bias_round %>% select(scenario_id, bias_chr),
                     relbias_round %>% select(relbias_chr),
                     empse_round %>% select(empse_chr),
                     modse_round %>% select(modse_chr))
  test3 <- plot_df3 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id)) 
  #View(test3)
  plot3 <- plot_df3 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id)) %>% 
    select(-(scenario_id)) %>%
    mutate(fontface_chr = "plain") %>%
    add_row(bias_chr = "Bias \n (MCSE)",
            relbias_chr = "Relative bias (%) \n (MCSE)",
            empse_chr = "Empirical SE \n (MCSE)",
            modse_chr = "Model SE \n (MCSE)",
            y_height = 0.5,
            fontface_chr = "bold") %>%
    ggplot(aes(x = 0, y = scenario_id)) +
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.y=element_blank() )+
    geom_text(aes(x= 0, y = y_height, label =  bias_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.075, y = y_height, label =  relbias_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.150, y = y_height, label =  empse_chr, fontface = fontface_chr), size = 2.5)+
    geom_text(aes(x= 0.215, y = y_height, label =  modse_chr, fontface = fontface_chr), size = 2.5)+
    scale_x_continuous("", limits=  c(-0.05, .255), expand = c(0, 0))+
    ylim(c(y_min,y_max))
  
  plot3
  
  #########################################
  # Plot 4. Forest/lollipop with Coverage
  #########################################
  
  
  plot_df4 <- all_res %>%
    filter(stat == "cover") %>%
    mutate(est_high = est + 1*mcse,
           est_low = est - 1*mcse) %>%
    select(est, mcse, est_low, est_high, scenario_id) %>%
    mutate(est_true = 0.95)
  
  test4 <- plot_df4 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id))
  #test4
  plot4 <- plot_df4 %>%
    filter(scenario_id %in% scen_df$scenario_id) %>%
    left_join(scen_df %>% select(scenario_id, y_height), 
              by = join_by(scenario_id)) %>%
    ggplot(aes(x = est, y = y_height)) +
    theme_classic()+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text( size = 6, angle = 45, vjust = 0.5, hjust=0.5),
          axis.title.x = element_text(face="bold", size = 8)) +
    geom_segment(x = 0.95, xend = 0.95, y = y_min-0.5, yend = y_max-1, colour = "#830051")+
    geom_segment(aes(x= est_true, xend = est, y = y_height, yend = y_height ),
                 alpha = 0.4, colour = "#830051")+
    # geom_errorbar(aes(xmin=est_low, xmax=est_high),
    #                width=1,                    # Width of the error bars
    #                position=position_dodge(.9),
    #                colour = "black",
    #                alpha = 0.6)+
    geom_point(aes(x=est_low, y = y_height), shape = 91, size = 2.5, colour = "#830051")+
    geom_point(aes(x=est_high, y = y_height), shape = 93, size = 2.5, colour = "#830051")+
    geom_point(alpha = 0.7, colour= "#830051") +
    scale_x_continuous("Coverage", limits = c(0.75, 1.00), labels = scales::percent,
                       breaks = seq(from = 0.75, to = 1.00, by = 0.05),
                       expand = c(0, 0))+
    ylim(c(y_min,y_max))
  plot4
  
  ##################################################
  ##################################################
  
  wid1 <- 2.6
  wid1null <- -0.4
  wid2 <- 2.0
  wid2null <- -0.2
  wid3 <- 6
  wid3null <- -0.2
  wid4 <- 2.0
  wid4null <- 0.2
  
  figure1 <- plot_grid(plot1, NULL,plot2, NULL, plot3, NULL, plot4, NULL, nrow = 1,
                       rel_widths = c(wid1, wid1null, wid2, wid2null, wid3,
                                      wid3null, wid4, wid4null),
                       align = "hv")
  
  figure1
  
  
  #################################
  ## getwd()
  if(!dir.exists("plots")){
    dir.create("plots")
  }  
  
  saveRDS(figure1, "plots/supfigure_bench_perform.RDS")
  
}


#####################################


figure_cetux <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", cetux_jobname, "/plots/supfigure_bench_perform.RDS"))
figure_niv <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", niv_jobname, "/plots/supfigure_bench_perform.RDS"))

figure_all <- plot_grid(NULL, figure_cetux, NULL, figure_niv,
                        rel_heights = c(0.03, 0.5,0.03,  0.5), 
                        labels = c("(a) Cetuximab OS", "", "(b) Nivolumab PFS", ""),
                        label_size = 10,
                        label_x = -0.062,
                        ncol = 1)
figure_all
#plot(0,0)

pdf(file = paste0(store_directory, "figures/figure_bench_perform.pdf"),   
    width = 8, 
    height = 6.3) 
figure_all

dev.off()

tiff(file = paste0(store_directory, "figures/figure_bench_perform.tiff"),   
     width = 8, 
     height = 6.3,
     units = 'in',  
     res = 1200, 
     compression = "lzw")

figure_all
dev.off()



############################################

scenarios_rstpm

res_test <- res 

res_test %>%
  filter(method  == "rstpm2",
         df == 6,
         isim != 0,
         value > 4)

res_test %>%
  filter(method  == "rstpm2",
         df == 6,
         isim == 789)

res_test[9842,]

summary(res_test )
summary(as.factor(res_test$method))
mean(res_test$value_se)
sd(res_test$value[res_test$value < 4])
sd(res_test$value[res_test$value < 6])
summary(res_test$value)

sum(res_test$value > 3.4)
res_test$value[res_test$value > 3.5]

##############################################

scenarios_rstpm
res_test2 <- readRDS(paste0("rstpm2_scen", 4, "_res.rds"))
res_test2 <- res_test2 %>%
    filter(isim == 819)
res_test2
res_test2 <- res_test2 %>%
  filter(isim == 785)
res_test2
res_test2 <- res_test2 %>%
  filter(isim == 789)
res_test2

res_test4 <- readRDS(paste0("rstpm2_scen", 4, "_est.rds"))
res_test4 %>%
  filter(isim == 785)

##############################################

res_test3 <- readRDS("scen_rstpm2_4/rstpm2_est785.rds")
res_test3



##############################


res_rstpm2 %>%
  filter(isim == 785)

res <- rbind(res_flex,
             res_rstpm2) %>%
  group_by(df, method) %>%
  mutate(scenario_id = cur_group_id())





names(sim_data1)
test_sim <- sim_data1 %>%
  filter(i == 785)
test_sim  


stpm2mod0 <- try(stpm2(Surv(time, event)~1, 
                       data = test_sim , 
                       df=4))


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

