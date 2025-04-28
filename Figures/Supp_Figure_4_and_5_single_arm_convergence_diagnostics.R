#####################
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
library(tibble)
library(tidyr)


jobname <- "cetux16"
#jobname <- "niv8" 


cetux_jobname <-  "cetux16"
niv_jobname <- "niv8" 

user <- Sys.info()["user"]
store_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim1_simulations/slurm/")

source("R/fit_model.R")
source("R/estimands.R")
source("R/visualise.R")
source("R/performance.R")



##############################################
setwd(paste0(store_directory , "simsurvextrap_slurm_", jobname, "/"))



scenarios_nsim <- read_csv(paste0("scenarios_nsim_", jobname, ".csv" ))
scenarios <- read_csv(paste0("scenarios_", jobname, ".csv" ))



###################################################
# Combine data together.
###################################################

warning_vec <- c("prop_divergent",
                  "warning_any",
                  "divergent",
                  "treedepth",
                  "BayesianFrac",
                  "SamplingProb",
                  "rhatHigh",
                  "BulkEssLow",
                  "TailEssLow")


for(i in scenarios$scenario_id){
 #i <- 1
  temp <- readRDS(paste0("scen", i, "_est.rds"))
  
  temp <- temp %>% 
    filter(estimand %in% warning_vec) %>%
    mutate(scenario_id = i) # %>%
  
  
  if(i == scenarios$scenario_id[1]){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }

  res <- data.frame(res)
  print(paste0("scenario ", i, "/", max(scenarios$scenario_id)))
}

dgm_num <- 1

scenarios_plot <- scenarios %>%
  filter(dgm_id == dgm_num) %>%
  mutate(stan_fit_method = if_else(stan_fit_method == "mcmc", "MCMC", "Opt"))

#names(scenarios_plot)
label_lookup <- scenarios_plot %>%
  select(c(scenario_id, mspline_df, prior_hsd_shape, 
           prior_hsd_rate,prior_hscale_mean, prior_hscale_sd ,stan_fit_method, 
           smooth_model, mspline_bsmooth, hscale_default)) %>%
  mutate(label_df = paste0("df~`=`~", mspline_df),
         label_gamma = paste0("sigma", "~`~`~", "`G`*`amma`*`(`*", 
                              prior_hsd_shape, "*`,`*" , 
                              prior_hsd_rate, "*`)`"),
         label_eta = paste0("eta", "~`~`~", "`N`*`(`*", 
                            prior_hscale_mean, "*`,`*" , 
                            prior_hscale_sd, "*`)`"),
         label_method = stan_fit_method) %>%
  mutate(label_df = fct_inorder(label_df),
         label_gamma = fct_inorder(label_gamma),
         label_eta = fct_inorder(label_eta),
         label_method = fct_inorder(label_method))


warning_res <- res %>%
  mutate(value = value*100) %>%
  left_join(scenarios_nsim %>% 
              mutate(num = row_number()) %>% 
              select(scenario_id, isim), 
            by = join_by(scenario_id, isim)) %>%
  group_by(scenario_id, estimand) %>%
  mutate(nsim = n()) %>%
  summarise(mean_value = mean(value),
            count = n()) %>%
  ungroup()  %>%
  left_join(label_lookup, by = join_by(scenario_id)) %>%
  filter(stan_fit_method == "MCMC") %>%
  mutate(warning_label = factor(estimand, 
                                levels = c("warning_any",
                                           "divergent",
                                           "rhatHigh",
                                           "BulkEssLow",
                                           "TailEssLow"),
                                labels = c(paste("Any warning messages", "reported", 
                                                 sep =  "\n"),
                                           paste("Divergent transitions", "reported", 
                                                 sep =  "\n"),
                                           "R-hat high",
                                           "bulk-ESS low",
                                           "tail-ESS low"))) %>%
  filter(!is.na(warning_label)) %>%
  mutate(x_null = 1)
#View(warning_res)


smooth_model_type  <-  "random_walk"

if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"

summary(as.factor(warning_res$estimand))
summary(as.factor(warning_res$warning_label))

warning_plot <- warning_res %>%
  filter(mspline_df %in% c(3,6,10),
         prior_hsd_rate %in% c(1,5,20),
         smooth_model %in% smooth_model_type,
         mspline_bsmooth %in% TRUE,
         hscale_default %in% TRUE) %>%
  mutate(label_df = fct_inorder(label_df),
         label_gamma = fct_inorder(label_gamma),
         label_method = fct_inorder(label_method)) %>%
  ggplot(aes(x = x_null, y = mean_value/100, fill = warning_label )) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y=element_text(size = 6),
        legend.position.inside = c(0.84, 0.89),
        strip.text = element_text(size=6),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8),
        axis.title=element_text(size=8))+
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Percentage with Stan warning message") +
  scale_fill_discrete(name = "Warning message")+
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  facet_wrap(~label_df+label_gamma, ncol = 3, labeller = label_parsed)

warning_plot
saveRDS(warning_plot, paste0("plots/sup_figure_warning", smooth_model_ext, ".rds"))

###################################
getwd()

smooth_model_type  <- "exchangeable"

if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"



figure1 <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", cetux_jobname, "/", 
                          "plots/sup_figure_warning", smooth_model_ext, ".rds"))
figure2 <- readRDS(paste0(store_directory , "simsurvextrap_slurm_", niv_jobname, "/",
                           "plots/sup_figure_warning", smooth_model_ext, ".rds"))


plot_all <- plot_grid(NULL,figure1, NULL, figure2, rel_heights = c(0.08, 0.5, 0.08, 0.5),  
                      labels = c("(a) Cetuximab-OS", "", "(b) Nivolumab-PFS", ""),
                      label_size = 10,
                      label_x = -0.092, ncol = 1,align = "v")

plot_all


#plot_all
legend1 <- cowplot::get_legend(figure1)
leg_gg <- as_ggplot(legend1)
leg_gg
#legend1


figure1_noleg <- figure1 + theme(legend.position='none')
figure2_noleg <- figure2 + theme(legend.position='none')

plot_warnings <- ggdraw(plot_grid(plot_grid(NULL,figure1_noleg, NULL, figure2_noleg, rel_heights = c(0.04, 0.5, 0.04, 0.5),  
                           labels = c("(a) Cetuximab OS", "", "(b) Nivolumab PFS", ""),
                           label_size = 10,
                           label_x = -0.092, ncol = 1,align = "v"),
                 plot_grid(NULL, plot_grid(legend1, NULL, ncol = 1, rel_heights = c(1, 100)), ncol=1),
                 rel_widths=c(1, 0.35)))
plot_warnings

jobname

pdf(file = paste0(store_directory, "figures/figure_warnings_", smooth_model_ext, ".pdf"),   
    width = 5.5, 
    height = 7) 

plot_warnings

dev.off()



tiff(file = paste0(store_directory, "figures/figure_warnings_", smooth_model_ext, ".tiff"),   
     width = 5.5, 
     height = 7,
     units = 'in',  
     res = 1200, 
     compression = "lzw")

plot_warnings

dev.off()


smooth_model_type

# extract results for text
View(warning_res %>%
  filter(mspline_df %in% c(10),
         prior_hsd_rate %in% c(1),
         smooth_model %in% "random_walk",
         mspline_bsmooth %in% TRUE,
         hscale_default %in% TRUE))


View(scenarios)
test <- readRDS(paste0("scen", 6, "_est.rds"))
summary(as.factor(test$estimand))
#diverge_count <- res
#summary(as.factor(diverge_count$estimand))
diverge_test <- test %>%
  filter(estimand == "divergent") %>%
  summarise(mean = 100*sum(value)/n())
diverge_test

sim_id <- test  %>%
  filter(estimand == "divergent", value > 0) %>%
  pull(isim)
#  filter(estimand == "prop_divergent") %>%

    
test_prop <- test %>%
    filter(estimand == "prop_divergent") %>%
    filter(isim %in% sim_id) %>%
  summarise(mean = 100*sum(value)/n()) %>%
  pull()
test_prop  
#mean(test_prop$value)*100
#summary(as.factor())

