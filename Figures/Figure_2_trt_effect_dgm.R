
#############################################################
# Plot data-generating mechanism survival/hazard functions.
#############################################################

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
library(emg)


# Folder where figures will be stored.
#store_directory <- "/home/lshit9/survival_extrapolation/simulation1/figures"
store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/simsurvextrap_trt3"
store_res <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/figures"


source("Functions/functions_treatment_effect.R")

# Import Bonner trial Cetuximab arm.
surv_df <- readRDS("Data/cetuximab_OS.rds")
# End of trial follow-up time for rmst computation
maxT_data <- max(surv_df$years)
k_true <-  3
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)

label_vec <- c("(a) Scenario 1: Constant effect (proportional hazards)",
               "(b) Scenario 2: Waning effect", 
               "(c) Scenario 3: Delayed then waning effect",
		"(d) Scenario 4: Crossing survival curves")

for(scenario_num in 1:4){
#scenario_num <- 1

big_df <-  readRDS(paste0(store_directory , "/", "dgm", scenario_num , "/sim_data_big.rds"))

km_fit <- survfit(Surv(time, event) ~ trt, data=big_df %>% slice(1:100000))

km_plot <- ggsurvplot(km_fit, data=big_df %>% slice(1:100000))

surv_df <- km_plot$data.survplot %>%
  mutate(trt = factor(trt, labels = c("Control", "Active"))) %>%
  arrange(trt) %>%
  mutate(random = runif(n(),0,1)) %>%
  filter(time <= 5) %>%
  arrange(random) %>%
  slice(1:100000)

plot_surv <- surv_df %>%
ggplot(aes(x = time, y = surv, colour = trt))+
  theme_classic()+
  theme(legend.position = c(0.84, 0.96),
        legend.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0, 0, 0, 0.1), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_step()+
  ylab("Overall survival (%)")+
  scale_y_continuous(lim = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)", lim = c(0, 5), breaks = 0:5)+
  scale_colour_discrete(name = NULL)

# log hazard ratio coefficients.
beta1 <- log(0.7)
beta2 <- c(log(0.5)/(1-tanh(-1.2)), 0.8, -1.2)
beta3 <- c(0, -2.8, 0.8, 0.4, 0.35)
beta4 <- c(5, 6, log(0.5), log(1.8))

hr <- scenario_hr(scenario_num)
haz <- scenario_haz(scenario_num)
beta_true <- get(paste0("beta", scenario_num))

hr_df <- tibble(time = seq(from = 0, to = 5, length.out = 1e4)) %>%
  mutate(hr = hr(t = time,
                 gamma = true_mod$coefficients,
                 knots = true_mod$knots, 
                 beta =  beta_true) )

plot_hr <- hr_df %>%
  ggplot(aes(x = time, y = hr))+
  theme_classic()+
  theme(plot.margin = unit(c(0, 0.2, 0, 0.15), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_hline(yintercept = 1, colour = "gray55", linetype = "dashed")+
  geom_line(col = "#830051")+
  ylab("Hazard ratio")+
  scale_y_log10(lim = c(0.42,2.14), 
                breaks = c(0.50, 1.00, 2.00))+
  scale_x_continuous("Time (years)", lim = c(0, 5), 
                     breaks = 0:5)+
  scale_colour_discrete(name = NULL)

figure <- plot_grid(plot_surv+
                      theme(plot.title = element_text(vjust = 5,
                                                      size=10, 
                                                      face="bold"),
                            plot.title.position = "plot")+
                      labs(title = label_vec[scenario_num]),
                    plot_hr, nrow = 1,align = "h")
  

assign(paste0("figure", scenario_num), figure)
saveRDS(figure, paste0(store_directory , "/plots/figure4_",scenario_num,".RDS"))
print(scenario_num)

}

plot_all <- plot_grid(NULL,figure1, NULL, figure2, NULL, figure3, NULL, figure4,
                      rel_heights = c(0.08, 0.5, 0.08, 0.5, 0.08, 0.5, 0.08, 0.5),  
                      ncol = 1,align = "v")


plot_all


pdf(file = paste0(store_res, "/figure_dgm_two_arm.pdf"),   
    width = 5.2, 
    height = 7.0) 
plot_all
dev.off()


tiff(file = paste0(store_res, "/figure_dgm_two_arm.tiff"),   
     width = 5.2, 
     height = 7.0,
     units = 'in',  
     res = 1200, 
     compression = "lzw")
plot_all
dev.off()

