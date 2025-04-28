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

######################################

store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")


# Extract Bonner trial Cetuximab arm.
surv_df <- cetux[cetux$treat=="Cetuximab",]
# End of trial follow-up time for rmst computation
maxT_data <- max(surv_df$years)
k_true <-  3
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
summary(true_mod,t = c(2,3,5), type = "survival")


km_fit <- survfit(Surv(years, d) ~ 1, data=surv_df)
km_plot <- ggsurvplot(km_fit,data=surv_df)

pred_haz <- standsurv(true_mod, 
                      t=seq(from = 0.01, to = 5, by=0.01),
                      type="hazard") 

pred_surv <- standsurv(true_mod, 
                       t=seq(from = 0.01, to = 5, by=0.01),
                       type="surv") 

plot_haz1 <- pred_haz %>%
  ggplot(aes(x = time, y = at1))+
  theme_classic()+
  theme(legend.position = c(0.84, 0.76),
        legend.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0, 0.2, 0, 0.15), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_line(linewidth  = 0.7)+
  scale_y_continuous("Hazard", lim = c(0,0.5)) +
  scale_x_continuous("Time (years)", lim = c(0, 5), breaks = 0:5) 

plot_surv1 <- pred_surv %>%
  ggplot(aes(x = time, y = at1))+
  theme_classic()+
  theme(legend.position = c(0.84, 0.76),
        legend.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0, 0, 0, 0.2), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_step(aes(x = time, y = surv), 
            km_plot$data.survplot, 
            colour = "deeppink3", 
            linewidth = 0.7,
            alpha = 0.9)+
  geom_line(linewidth  = 0.7, alpha = 0.6)+
  ylab("Overall survival (%)")+
  scale_y_continuous(lim = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)", lim = c(0, 5), breaks = 0:5)+
  theme(plot.title = element_text(vjust = 3.5,
                                  size=10, 
                                  face="bold"),
        plot.title.position = "plot")+
  labs(title = "(a) Cetuximab OS")

plot_grid(plot_surv1, plot_haz1, ncol = 1,align = "v")


# Nivolumav pfs
surv_df2 <- read_csv(here::here("nivo_data/nivolumab_digiKM_pfs.csv"),
                    show_col_types = FALSE)  %>%
  filter(treat==2) %>%
  mutate(years = time/12) %>%
  rename(d = status)

# end of trial follow-up time for rmst computation

maxT_data2 <- max(surv_df2$years)

k_true2 <-  6
true_mod2 <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df2, k=k_true2)

summary(true_mod2,t = 5, type = "rmst")
summary(true_mod2,t = 5, type = "survival")


km_fit2 <- survfit(Surv(years, d) ~ 1, data=surv_df2)
km_plot2 <- ggsurvplot(km_fit2,data=surv_df2)

pred_haz2 <- standsurv(true_mod2, 
                      t=c(seq(from = 0.01, to = 0.5, length.out = 1e5),
                          seq(from = 0.5, to = 5, length.out = 1e4)
                      ),
                      type="hazard") 

pred_surv2 <- standsurv(true_mod2, 
                       t=seq(from = 0.01, to = 5, length.out = 1e4),
                       type="surv") 

plot_haz2 <- pred_haz2 %>%
  ggplot(aes(x = time, y = at1))+
  theme_classic()+
  theme(legend.position = c(0.84, 0.76),
        legend.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0, 0.2, 0, 0.15), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_line(linewidth  = 0.7)+
  scale_y_continuous("Hazard", lim = c(0,3)) +
  scale_x_continuous("Time (years)", lim = c(0, 5), breaks = 0:5) 

plot_surv2 <- pred_surv2 %>%
  ggplot(aes(x = time, y = at1))+
  theme_classic()+
  theme(legend.position = c(0.84, 0.76),
        legend.title = element_blank(), 
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"),
        plot.margin = unit(c(0, 0, 0, 0.2), "cm"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=9))+
  geom_step(aes(x = time, y = surv), 
            km_plot2$data.survplot, 
            colour = "deeppink3", 
            linewidth = 0.7,
            alpha = 0.9)+
  geom_line(linewidth  = 0.7, alpha = 0.6)+
  scale_y_continuous("Progression free survival (%)", lim = c(0,1), labels = scales::percent)+
  scale_x_continuous("Time (years)", lim = c(0, 5), breaks = 0:5)+
  theme(plot.title = element_text(vjust = 3.5,
                                  size=10, 
                                  face="bold"),
        plot.title.position = "plot")+
  labs(title = "(b) Nivolumab PFS")


plot1 <- plot_grid(plot_surv1,plot_haz1, nrow = 1,align = "h")
plot2 <- plot_grid(plot_surv2,plot_haz2, nrow = 1,align = "h")
plot1
plot2

plot_all <- plot_grid(NULL,plot1, NULL, plot2, rel_heights = c(0.08, 0.5, 0.08, 0.5),  
          labels = c("(a) Cetuximab OS", "", "(b) Nivolumab PFS", ""),
          label_size = 10,
          label_x = -0.092, ncol = 1,align = "v")

plot_all <- plot_grid(NULL,plot1, NULL, plot2, rel_heights = c(0.05, 0.5, 0.05, 0.5),  
                      ncol = 1,align = "v")

plot_all

pdf(file = paste0(store_directory, "figures/figure_dgm_single_arm.pdf"),   
    width = 6, 
    height = 4.9) 
plot_all
dev.off()


tiff(file = paste0(store_directory, "figures/figure_dgm_single_arm.tiff"),   
     width = 6, 
     height = 4.9,
     units = 'in',  
     res = 1200, 
     compression = "lzw")
plot_all
dev.off()

##############################################
# number of events.


jobname <- "master_cetux10"

head(sim_data1)
summary(sim_data1$id)
summary(sim_data1$i)

sim_data1 %>%
  group_by(i) %>%
  summarise(count = n(),
            events = sum(event)) %>%
  summarise(mean = mean(events),
            sd = sd(events))# %>%
  pull()
