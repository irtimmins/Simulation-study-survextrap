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


jobname_cetux <- "cetux_OS"
jobname_niv <-  "niv_PFS" 
store_directory_base <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"
store_directory_cetux <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/simsurvextrap_cetux_short1/"
store_directory_niv <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/simsurvextrap_niv_short1/"

source("Functions/estimands.R")
source("Functions/performance.R")

##########################################
# Plot RMST at 5-y performance measures
##########################################

cetux_jobname <-  "cetux_OS"
niv_jobname <- "niv_PFS" 

est_id_choose <- "rmst1"
est_id_years <- 2

# also run for 
est_id_choose <- "rmst2"
est_id_years <- 3

for(job in c("cetux", "niv")){
  
	#job <- "cetux"
  if(job == "cetux"){
    store_directory <- store_directory_cetux 
    jobname <- jobname_cetux
	} else if (job == "niv") {
    store_directory <- store_directory_niv
    jobname <- jobname_niv

    }


  #jobname <- cetux_jobname
  setwd(store_directory)
  print(getwd())

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
    
    temp <- readRDS(paste0("scen", i, "_res.rds"))
    
    temp <- temp %>% 
      filter(estimand == "rmst") %>%
      mutate(scenario_id = i) # %>%
    
    
    if(i == scenarios$scenario_id[1]){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
    res <- data.frame(res)
    res <- res %>%
      filter(!is.na(value)) %>%
      filter(!is.na(value_se) | isim == 0)
  }
  saveRDS(res, "plots/res.RDS")
  
  
  perform_res <- NULL
  for(i in 1:nscen){
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
  saveRDS(perform_res, "plots/perform_res.RDS")
  
  all_res <- 
    perform_res %>%
    filter(stat %in% c( "nsim", "bias", "rbias","empse", "mse", "modelse", "cover")) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round))
  
  saveRDS(all_res, "plots/all_res.RDS")
  
  # convert relative bias to percentage
  # all_res$est_round[all_res$stat == "rbias"] <- 100*all_res$est[all_res$stat == "rbias"]
  # all_res$mcse[all_res$stat == "rbias"] <- 100*all_res$mcse[all_res$stat == "rbias"]
  # all_res$est_round[all_res$stat == "rbias"] <- sprintf("%.1f", round(all_res$est[all_res$stat == "rbias"] , 3))
  # all_res$mcse_round[all_res$stat == "rbias"] <- sprintf("%.2f", round(all_res$mcse[all_res$stat == "rbias"] , 3))
  # all_res$est_round[all_res$stat == "rbias"] <- paste0(all_res$est_round[all_res$stat == "rbias"], "%")
  # all_res$mcse_round[all_res$stat == "rbias"] <- paste0(all_res$mcse_round[all_res$stat == "rbias"], "%")
  
}

####################################
# Plot. 
####################################

#names(scenarios)
#summary(as.factor(scenarios$prior_hsd_rate))

for(job in c("cetux", "niv")){
  
  #job <- "cetux"


  if(job == "cetux"){
    store_directory <- store_directory_cetux 
    jobname <- jobname_cetux
	} else if (job == "niv") {
    store_directory <- store_directory_niv
    jobname <- jobname_niv

    }

  setwd(store_directory)
  print(getwd())

  for(stan_fit in c("mcmc", "opt")){
   
	

    #stan_fit <- "mcmc"
    
    #dir.create("plots")
    
    scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
    nscen <- nrow(scenarios)
    
    res <- readRDS("plots/res.RDS")
    perform_res <- readRDS("plots/perform_res.RDS")
    all_res <- readRDS("plots/all_res.RDS")
    
    for(smooth_model_type in c("random_walk", "exchangeable")){
      
      # smooth_model_type <- "random_walk"
      
      default_id <- scenarios %>%
        filter(mspline_df == 10,
               prior_hsd_rate == 1,
               stan_fit_method == "mcmc",
               smooth_model == "random_walk",
               mspline_bsmooth == T,
		maxT == 2) %>%
        pull(scenario_id)
      
      scen_df <- scenarios %>%
        filter(mspline_df %in% c(3,6, 10),
               prior_hsd_rate %in% c(1,5,20),
               stan_fit_method %in% stan_fit,
               smooth_model %in% smooth_model_type,
               mspline_bsmooth %in% TRUE,
               hscale_default %in% TRUE,
		maxT == 2) %>%
        mutate(label_hsd = paste0("`G`*`amma`*`(`*", 
                                  prior_hsd_shape, "*`,`*" , 
                                  prior_hsd_rate, "*`)`")) %>%
        select(scenario_id,
               mspline_df,
               label_hsd) %>%
        mutate(y_height = -row_number()) %>%
        mutate(mspline_df = as.character(mspline_df)) %>%
        mutate(fontface_chr = "plain") %>%
        add_row(tibble_row(scenario_id = 0, 
                           mspline_df = "M-spline df    ", 
                           label_hsd = paste0("bold(sigma", "~", "prior)"), 
                           y_height = 0,
                           fontface_chr = "bold")) %>%
        mutate(y_height = y_height + 10) 
      
      
      y_limits <- c(min(scen_df %>% pull(y_height)),
                    max(scen_df %>% pull(y_height))+0.1)

      #########################################
      # Plot 1. Model
      #########################################
      
      plot1 <- scen_df %>%
        mutate(default_label = if_else(scenario_id == default_id, "#830051", NA)) %>%
        ggplot(aes(x = 0, y = y_height)) +
        theme_classic()+
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks=element_blank(),
              axis.title.y=element_blank(),
              plot.margin = unit(c(0.1, 0, 0.2, 0.2), "cm"))+
        geom_text(aes(x= 0, y = y_height, label =  mspline_df, fontface = fontface_chr), size = 2.5)+
        geom_text(aes(x= 0.075, y = y_height, label =  label_hsd, fontface = fontface_chr), parse = T,size = 2.5)+
        #geom_label(aes(x= 0, y = y_height, label =  mspline_df, fontface = fontface_chr), size = 2.5, 
        #          data = subset(scen_df,scenario_id == default_id),
        #          fill = "#830051", alpha = 0.15)+
        #geom_label(aes(x= 0.075, y = y_height, label =  label_hsd, fontface = fontface_chr),
        #          fill = "#830051",parse = T, size = 2.5, alpha = 0.15, 
        #           data = subset(scen_df,scenario_id == default_id))+
        geom_label(aes(x= 0.0517 , y = y_height, 
                       label =  "                                 ", 
                       fontface = fontface_chr), size = 2.5, 
                   data = subset(scen_df,scenario_id == default_id),
                   colour = "#830051",
                   fill = "#830051", alpha = 0.15)+ 
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
                        bias %>% select(-c(scenario_id)))
      
      #View(test)
      if(est_id_choose == "rmst1"){
        if(job == "cetux") {
          xlim <- c(1.46, 1.80)
          xbreaks <- seq(from = 1.50, to = 1.76, by = 0.05)             
        } else{
          xlim <- c(0.86, 1.20) 
          xbreaks <- seq(from = 0.90, to = 1.16, by = 0.05)    
        }
      } else if(est_id_choose == "rmst2"){
        
        if(job == "cetux") {
          xlim <- c(2.01, 2.35)
          xbreaks <- seq(from = 2.05, to = 2.31, by = 0.05)             
        } else{
          xlim <- c(1.21, 1.55) 
          xbreaks <- seq(from = 1.25, to = 1.51, by = 0.05)    
        }
        
        
      }
        
      
      print(plot_df2 %>%
        filter(scenario_id %in% scen_df$scenario_id) %>%
        left_join(scen_df %>% select(scenario_id, y_height), 
                  by = join_by(scenario_id)) %>%
        pull(estimand_mean))
      
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
              axis.title.x = element_text( face="bold",size = 8),
              plot.margin = unit(c(0.1, 0, 0.2, 0), "cm")) + 
        geom_segment(x = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                     xend = plot_df2 %>% select(estimand_true) %>% slice(1) %>% pull(), 
                     y = 0.5, 
                     yend = 9.5, colour = "#830051")+
        geom_segment(aes(x= estimand_true, xend = estimand_mean, y = y_height, yend = y_height),
                     alpha = 0.4, colour = "#830051")+
        geom_point(aes(x=estimand_true+est_low, y = y_height), shape = 91, size = 2.5, colour = "#830051")+
        geom_point(aes(x=estimand_true+est_high, y = y_height), shape = 93, size = 2.5, colour = "#830051")+
        geom_point(alpha = 0.7, colour= "#830051") +
        scale_x_continuous(paste0("RMST at ",est_id_years,"-y"), limits = xlim,
                           breaks = xbreaks,
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
                modse_chr = "Model SD \n (MCSE)",
                y_height = 10,
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
              plot.margin = unit(c(0.1, 0, 0.2, 0), "cm"))+
        geom_text(aes(x= 0, y = y_height, label =  bias_chr, fontface = fontface_chr), size = 2.5)+
        geom_text(aes(x= 0.075, y = y_height, label =  relbias_chr, fontface = fontface_chr), size = 2.5)+
        geom_text(aes(x= 0.150, y = y_height, label =  empse_chr, fontface = fontface_chr), size = 2.5)+
        geom_text(aes(x= 0.215, y = y_height, label =  modse_chr, fontface = fontface_chr), size = 2.5)+
        scale_x_continuous("", limits=  c(-0.05, .255), expand = c(0, 0))+
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
        mutate(est_true = 0.95)
      
      test4 <- plot_df4 %>%
        filter(scenario_id %in% scen_df$scenario_id) %>%
        left_join(scen_df %>% select(scenario_id, y_height), 
                  by = join_by(scenario_id))
      test4
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
              axis.title.x = element_text(face="bold", size = 8),
              plot.margin = unit(c(0.1, 0.2, 0.2, 0), "cm")) +
        geom_segment(x = 0.95, xend = 0.95, y = 0.5, yend = 9.5, colour = "#830051")+
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
        ylim(y_limits)
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
      
      if(est_id_choose == "rmst1"){
        label_vec <- c("(a)", "(c)")
      } else(
        label_vec <- c("(b)", "(d)")
      )
      
      
      if(jobname == cetux_jobname) {
        figure_label <- paste0(label_vec[1]," Cetuximab OS, ", est_id_years, "-year data cut")
      } else{
        figure_label <- paste0(label_vec[2]," Nivolumab PFS, ", est_id_years, "-year data cut")
      }
      
      
      figure_title <- ggdraw() + 
        draw_label(
          figure_label,
          hjust = 0,
          y = 0.4,
          x = 0,
          size=10, 
          fontface="bold"
        ) +
        theme(
          plot.margin = margin(0, 0, 0, 7)
        )

      figure_row <- plot_grid(plot1, NULL,plot2, NULL, plot3, NULL, plot4, NULL, nrow = 1,
                           rel_widths = c(wid1, wid1null, wid2, wid2null, wid3,
                                          wid3null, wid4, wid4null),
                           align = "h")
      
      figure1   <- plot_grid(figure_title, figure_row,
                                      ncol = 1,
                                      rel_heights = c(0.08, 1))
      
      figure1
      #################################
      ## getwd()
      
      if(!dir.exists("plots")){
        dir.create("plots")
      }  
      
      
      if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
      if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
      
      print(getwd())
      #print(paste0("plots/supfigure_opt_perform_", smooth_model_ext,".RDS"))
      saveRDS(figure1, paste0(store_directory, "plots/figure_short_",stan_fit,"_perform_", smooth_model_ext,"_", est_id_choose,".RDS"))
    }
    
    
  }
}

#####################################

for(stan_fit in c("mcmc", "opt")){
  
  for(smooth_model_type in c("random_walk","exchangeable")){ 
    
    if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
    if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
    
    figure_cetux_rmst1 <- readRDS(paste0(store_directory_cetux , "plots/figure_short_", stan_fit, "_perform_",smooth_model_ext,"_rmst1.RDS"))
    figure_cetux_rmst2 <- readRDS(paste0(store_directory_cetux , "plots/figure_short_", stan_fit, "_perform_",smooth_model_ext,"_rmst2.RDS"))
    
    figure_niv_rmst1 <- readRDS(paste0(store_directory_niv , "plots/figure_short_", stan_fit, "_perform_",smooth_model_ext,"_rmst1.RDS"))
    figure_niv_rmst2 <- readRDS(paste0(store_directory_niv , "plots/figure_short_", stan_fit, "_perform_",smooth_model_ext,"_rmst2.RDS"))

    figure_all <- plot_grid(figure_cetux_rmst1,
                            figure_cetux_rmst2,
                            figure_niv_rmst1,
                            figure_niv_rmst2,
                            ncol = 1)
    
    figure_all  
    
    
    pdf(file = paste0(store_directory_base, "figures/sup_figure_short_", stan_fit,"_perform_",smooth_model_ext, ".pdf"),   
        width = 8, 
        height = 10.5) 
    print(figure_all)
    dev.off()
    
    tiff(file = paste0(store_directory_base, "figures/sup_figure_short_",stan_fit,"_perform_",smooth_model_ext,".tiff"),   
         width = 8, 
         height = 10.5,
         units = 'in',  
         res = 1200, 
         compression = "lzw")
    print(figure_all)
    dev.off()
    
    
    
  }
}





