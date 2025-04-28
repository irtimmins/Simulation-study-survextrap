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
# Combine results across iterations.
##########################################

#jobname <- "master_cetux12"
#jobname <- "master_niv4_smooth" 
#jobname <- "master_cetux10"  
#jobname <- "master_niv3"

jobname <-  "cetux16"
#jobname <- "niv8"


  setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))
  
  #dir.create("plots")
  
  scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))
  nscen <- nrow(scenarios)
  nsim <- 1000

for(i in scenarios$scenario_id){
  
  est <- NULL
  
  for(j in 1:nsim){
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", scenario ", i))
    if(file.exists(paste0("scen", i, "/est", j,".rds"))){
      temp <- readRDS(paste0("scen", i, "/est", j,".rds"))  %>%
        filter(estimand == "rmst") %>%
        mutate(isim = j, mod_type = "sim") 

      
      est <- rbind(est, temp)
    }
    
    
  }
  #nrow(est)
  saveRDS(est, paste0("scen", i, "_est_rmst.rds"))
  
  true <- readRDS("true_est.RDS")
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true") %>%
    filter(estimand == "rmst")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("scen", i, "_res_rmst.rds"))
  
}

#head(res)

##########################################
# Plot RMST at 5-y performance measures
##########################################

#cetux_vec <- c("master_cetux10","master_cetux12")
#niv_vec <- c("master_niv3","master_niv4_smooth")
cetux_vec <- "cetux16"
niv_vec <- "niv8" 
bsmooth_vec <- c(cetux_vec[1], niv_vec[1])  

#eval(paste0(scen_name, "_vec"))
scen_name <- "cetux"
job_vec <- get(paste0(scen_name, "_vec"))
for(jobname in  job_vec){
  #jobname <- niv_vec[1]
  #jobname <- "cetux16"
  #jobname <- "master_niv4_smooth"
  
  setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))
  
  #dir.create("plots")
  
  scenarios <- read_csv(paste0("scenarios_", jobname ,".csv"))

  scenarios <- scenarios %>%
    filter(prior_hscale_mean == 0)
  #nscen <- nrow(scenarios)
  
  #res_test <- readRDS("scen66_est_rmst.rds")
  #res_test
  #View(res_test)
  #summary(as.factor(res_test$estimand_id))
  
  ###################################################
  # Combine data together.
  ###################################################
  
  for(i in scenarios$scenario_id){
    
    temp <- readRDS(paste0("scen", i, "_res_rmst.rds"))
    
    temp <- temp %>% 
      mutate(scenario_id = i) # %>%
    
    
    if(i == scenarios$scenario_id[1]){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
    res <- data.frame(res)
    res <- res %>%
      filter(!is.na(value))
  }
  #scenarios <- read_csv(file = store_scenarios)
  temp_all_res <- res %>%
    left_join(scenarios %>% 
                select(scenario_id, mspline_bsmooth), 
              by = "scenario_id") %>%
    filter(t == 5) %>%
    mutate(estimand_id = "rmst1")
  
  #est_id_choose <- "rmst1"
  perform_res <- NULL
  for(i in scenarios$scenario_id){

    temp_res <- sim_perform(res %>%
                              filter(scenario_id == i) %>%
                              filter(t == 5) %>%
                              mutate(estimand_id = "rmst1"), 
                            estimands = c("rmst")) %>%
      mutate(scenario_id = i) 
     
    
    
    #head(temp_res)
    #names(scenarios)
    
    # if(jobname == "master_niv4_smooth") {
    #   temp_res <- temp_res %>%
    #     mutate(bsmooth = "F",
    #            new_id = paste0(scenario_id, "_", bsmooth))
    # } else {
    #   temp_res <- temp_res %>%
    #     mutate(bsmooth = "T",
    #            new_id = paste0(scenario_id, "_", bsmooth))
    #   
    # }
    # 
    if(i == 1){ 
      perform_res <- temp_res
    } else {
      perform_res <- rbind(perform_res, temp_res)
    }
    
  }
  #head(perform_res)
  perform_res <- perform_res %>%
    left_join(scenarios %>% 
              select(scenario_id, mspline_bsmooth), 
            by = "scenario_id")
  
  #head(perform_res)
  #summary(as.factor(perform_res$estimand_id))
  temp_all_perform_res <- 
    perform_res %>%
    filter(stat %in% c( "nsim", "bias", "rbias","empse", "mse", "modelse", "cover")) %>%
    mutate(est_round = sprintf("%.3f", round(est, 3)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round))
  #summary(as.factor(temp_all_perform_res$estimand_id))
  if(jobname %in% bsmooth_vec) {
    all_res <- temp_all_res
    all_perform_res <- temp_all_perform_res
    all_scenarios <- scenarios
  } else {
    all_res <- rbind(all_res, temp_all_res)
    all_perform_res <- rbind(all_perform_res, temp_all_perform_res)
    all_scenarios <- rbind(all_scenarios, scenarios)
  }
  
}

saveRDS(all_res, paste0("../one_arm_bsmooth/all_res_", scen_name,".RDS"))
saveRDS(all_perform_res, paste0("../one_arm_bsmooth/all_perform_res_", scen_name,".RDS"))
saveRDS(all_scenarios, paste0("../one_arm_bsmooth/all_scenarios_", scen_name ,".RDS"))

#  saveRDS(all_res, "../one_arm_bsmooth/all_res_niv.RDS")
#  saveRDS(all_perform_res, "../one_arm_bsmooth/all_perform_res_niv.RDS")
#  saveRDS(all_scenarios, "../one_arm_bsmooth/all_scenarios_niv.RDS")

#head(all_scenarios)

#head(all_res)
#head(all_res)
#View(all_res)

####################################
# Plot.
####################################


test_jobname <- "cetux16"
#test_jobname <- "niv8"
#setwd()
test1 <- read_csv(paste0(store_directory , "simsurvextrap_slurm_", test_jobname, "/",
                         "scenarios_", test_jobname ,".csv"))


summary(as.factor(test1$mspline_bsmooth))
View(test1)

setwd(paste0(store_directory , "one_arm_bsmooth", "/"))

smooth_model_type <- "random_walk"
#smooth_model_type <- "exchangeable"

if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"


for(scen_name in c("cetux", "niv")){
#scen_name <- "cetux"
all_scenarios <- readRDS(paste0("all_scenarios_", scen_name, ".RDS"))
all_perform_res <- readRDS(paste0("all_perform_res_", scen_name, ".RDS"))
all_res <- readRDS(paste0("all_res_", scen_name, ".RDS"))

#

all_scenarios <- all_scenarios %>%
  mutate(new_id = paste0(scenario_id,"_",  substr(mspline_bsmooth, start = 1, stop = 1)))

all_perform_res <- all_perform_res %>%
  mutate(new_id = paste0(scenario_id,"_",  substr(mspline_bsmooth, start = 1, stop = 1)))

all_res <- all_res %>%
  mutate(new_id = paste0(scenario_id,"_",  substr(mspline_bsmooth, start = 1, stop = 1)))

#length(unique(all_scenarios$new_id))
#length(all_scenarios$new_id)
#names(all_res)
#View(all_scenarios)
#summary(as.factor(all_scenarios$mspline_bsmooth))

scen_df <- all_scenarios  %>%
  filter(mspline_df %in% c(6, 10),
         prior_hsd_rate %in% c(1,5,20),
         prior_hscale_mean == 0,
         hscale_default %in% TRUE,
         stan_fit_method == "mcmc",
         smooth_model %in% smooth_model_type) %>%
  mutate(label_hsd = paste0("`G`*`amma`*`(`*", 
                            prior_hsd_shape, "*`,`*" , 
                            prior_hsd_rate, "*`)`"))  %>%
  arrange(!mspline_bsmooth) %>%
  arrange(prior_hsd_rate) %>%
  arrange(mspline_df) %>%
    select(new_id,
           scenario_id,
         mspline_df,
         label_hsd,
         mspline_bsmooth) %>%
  mutate(y_height = -row_number()) %>%
  mutate(mspline_df = as.character(mspline_df),
         mspline_bsmooth = as.character(mspline_bsmooth)) %>%
  mutate(mspline_bsmooth = str_to_title(mspline_bsmooth)) %>%
  mutate(mspline_bsmooth = if_else(mspline_bsmooth == "True", "Smoothed", "Standard")) %>%
  mutate(fontface_chr = "plain") %>%
  add_row(tibble_row(scenario_id = 0, 
                     mspline_df = "M-spline df", 
                     mspline_bsmooth = "Basis",
                     label_hsd = paste0("bold(sigma", "~", "prior)"), 
                     y_height = 0.5,
                     fontface_chr = "bold")) #%>%
#  mutate(y_height = y_height + 10) 

#scen_df$mspline_bsmooth <- if_else(scen_df$mspline_bsmooth == "True", "Smoothed", "Unsmoothed")
#head(scen_df)
#tail(scen_df)
#View(scen_df)

y_max <- max(scen_df %>%
               select(y_height) %>%
               pull())

y_min <- min(scen_df %>%
               select(y_height) %>%
               pull())


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
  geom_text(aes(x= 0, y = y_height, label =  mspline_bsmooth, fontface = fontface_chr), parse = F, size = 2.5)+
  geom_text(aes(x= 0.06, y = y_height, label =  mspline_df, fontface = fontface_chr), parse = F,size = 2.5)+
  geom_text(aes(x= 0.12, y = y_height, label =  label_hsd, fontface = fontface_chr), parse = T,size = 2.5)+
  scale_x_continuous("", limits = c(-0.05, .17), expand = c(0, 0))+
  ylim(c(y_min,y_max))

plot1

#########################################
# Plot 2. Forest/lollipop with RMST
#########################################

est_id_choose <- "rmst1"
#summary(as.factor(all_res$estimand_id))
est_mean <- all_res %>%
  filter(estimand_id == est_id_choose) %>%
  filter(isim != 0) %>%
  group_by(new_id) %>%
  summarise(estimand_mean = mean(value))

#head(all_res)
#head(est_mean)

est_true <- all_res %>%
  filter(estimand_id == est_id_choose) %>%
  filter(isim == 0) %>%
  group_by(new_id) %>%
  summarise(estimand_true = mean(value)) %>%
  ungroup() 

#head(est_true)
#head(bias)

bias <- all_perform_res %>%
  filter(stat == "bias") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id)


plot_df2 <- est_mean %>%
  left_join(est_true, join_by(new_id)) %>%
  left_join(bias, join_by(new_id))

#View(test)

if(scen_name== "cetux") {
  xlim <- c(3.10, 3.29) 
  xbreaks <- seq(from = 3.10, to = 3.29, by = 0.05)                      
} else{
  xlim <- c(1.90, 2.09)
  xbreaks <- seq(from = 1.90, to = 2.09, by = 0.05)
}

plot2 <- plot_df2 %>%
  filter(new_id %in% scen_df$new_id) %>%
  left_join(scen_df %>% select(new_id, y_height), 
            by = join_by(new_id)) %>%
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
               y = y_min-1, 
               yend = y_max-1, colour = "#830051")+
  geom_segment(aes(x= estimand_true, xend = estimand_mean, y = y_height, yend = y_height),
               alpha = 0.4, colour = "#830051")+
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

bias_round <- all_perform_res %>%
  filter(stat == "bias") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%                                                    
  mutate(bias_chr = paste0(est_round, " (", mcse_round, ")")) 

relbias_round <- all_perform_res %>%
  filter(stat == "rbias") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_round = sprintf("%.2f", round(100*est, 2)),
         mcse_round = sprintf("%.2f", round(100*mcse, 2)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%
  mutate(relbias_chr = paste0(est_round, "%", " (", mcse_round, "%)"))


#summary(as.factor(all_res$stat))
empse_round <- all_perform_res %>%
  filter(stat == "empse") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_round = sprintf("%.3f", round(est, 3)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%
  mutate(empse_chr = paste0(est_round, " (", mcse_round, ")"))

modse_round <- all_perform_res %>%
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

#head(plot_df3)  

# plot_df3 <- cbind( bias_round %>% select(scenario_id, bias_chr),
#                    relbias_round %>% select(relbias_chr),
#                    empse_round %>% select(empse_chr),
#                    modse_round %>% select(modse_chr))


# test3 <- plot_df3 %>%
#   filter(scenario_id %in% scen_df$scenario_id) %>%
#   left_join(scen_df %>% select(scenario_id, y_height), 
#             by = join_by(scenario_id)) 
#View(test3)

plot3 <- plot_df3 %>%
  filter(new_id %in% scen_df$new_id) %>%
  left_join(scen_df %>% select(new_id, y_height), 
            by = join_by(new_id)) %>% 
  select(-(new_id)) %>%
  mutate(fontface_chr = "plain") %>%
  add_row(bias_chr = "Bias \n (MCSE)",
          relbias_chr = "Relative bias (%) \n (MCSE)",
          empse_chr = "Empirical SE \n (MCSE)",
          modse_chr = "Model SD \n (MCSE)",
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
  scale_x_continuous("", limits=  c(-0.05, .265), expand = c(0, 0))+
  ylim(c(y_min,y_max))

plot3

#########################################
# Plot 4. Forest/lollipop with Coverage
#########################################

plot_df4 <- all_perform_res %>%
  filter(stat == "cover") %>%
  mutate(est_high = est + 1*mcse,
         est_low = est - 1*mcse) %>%
  select(est, mcse, est_low, est_high, new_id) %>%
  mutate(est_true = 0.95)

#head(all_perform_res)
# test4 <- plot_df4 %>%
#   filter(scenario_id %in% scen_df$scenario_id) %>%
#   left_join(scen_df %>% select(scenario_id, y_height), 
#             by = join_by(scenario_id))
#test4

plot4 <- plot_df4 %>%
  filter(new_id %in% scen_df$new_id) %>%
  left_join(scen_df %>% select(new_id, y_height), 
            by = join_by(new_id)) %>%
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
  geom_segment(x = 0.95, xend = 0.95, y = y_min-1, yend = y_max-1, colour = "#830051")+
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
  scale_x_continuous("Coverage", limits = c(0.8, 1.00), labels = scales::percent,
                     breaks = seq(from = 0.8, to = 1.00, by = 0.05),
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

figure1 <- plot_grid(plot1, NULL,plot2, NULL, plot3, NULL, plot4, NULL, nrow = 1,
                     rel_widths = c(wid1, wid1null, wid2, wid2null, wid3,
                                    wid3null, wid4, wid4null),
                     align = "h")

figure1


#################################
## getwd()
#scen_name

saveRDS(figure1,  paste0("figure_smooth_", scen_name, "_",smooth_model_ext,".RDS"))

#####

}


#####################################


figure_cetux <- readRDS(paste0("figure_smooth_", "cetux", "_",smooth_model_ext,".RDS"))
figure_niv <- readRDS(paste0("figure_smooth_", "niv", "_",smooth_model_ext,".RDS"))

figure_all <- plot_grid(NULL, figure_cetux, NULL, figure_niv,
                        rel_heights = c(0.03, 0.5,0.03,  0.5), 
                        labels = c("(a) Cetuximab OS", "", "(b) Nivolumab PFS", ""),
                        label_size = 10,
                        label_x = -0.062,
                        ncol = 1)
figure_all
#plot(0,0)

pdf(file = paste0(store_directory, "figures/figure_mcmc_bsmoothF_perform_",
                  smooth_model_ext,".pdf"),   
    width = 9, 
    height = 7.3) 
figure_all
dev.off()

tiff(file = paste0(store_directory, "figures/figure_mcmc_bsmoothF_perform_",
                   smooth_model_ext,".tiff"),   
     width = 9, 
     height = 7.3,
     units = 'in',  
     res = 1200, 
     compression = "lzw")

figure_all
dev.off()




