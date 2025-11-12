#############################
library(survextrap)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(cowplot)
library(purrr)
library(posterior)
library(readr)
library(purrr)
library(forcats)
library(ggtext)

store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/figures"

##########################################
# Test cases.
##########################################

surv_df <- readRDS("Data/cetuximab_OS_simulated_example.rds")
surv_df <- readRDS("Data/nivolumab_PFS.rds")


k <- 2
prior_choice <- c(1, 5, 20)[k]
surv_df %>% rename(event = d, time = years)
# mspline_evenly_spaced = list(knots = seq(from = 0.5, to = 5, length.out = 9))

nd_mod <- survextrap(Surv(time, event) ~ 1, 
                     data=surv_df,#%>% 
                       #rename(event = d, time = years), 
                     smooth_model =  "random_walk",
                     prior_hsd = p_gamma(2, prior_choice)) #, mspline = mspline_evenly_spaced)

#nd_mod$mspline$knots
#?p_gamma

prior_df <- tibble(x = seq(from = 0, to = 10, length.out = 100)) %>%
  mutate(dist = dgamma(x = x, shape = 2, rate= prior_choice)) %>%
  mutate(type = "Prior")

#prior_df
prior_df %>%
  ggplot(aes(x = x, y = dist))+
  theme_classic()+
  geom_line()

model_draws <- as_draws_df(nd_mod$stanfit) %>%
  as_tibble() %>%
  select(starts_with("hsd")) %>%
  rename(sigma = `hsd[1]`)

post_dens <- density(model_draws$sigma, bw = 0.15)
post_df <- tibble(x = post_dens$x,
                  dist = post_dens$y) %>%
  mutate(type = "Posterior")

# post_df
# post_dens$x
#plot(post_dens)
#names(model_draws)
#head(model_draws)

#summary(model_draws)
post_df  %>%
  bind_rows(prior_df) %>%
ggplot(aes(x=x, y = dist, colour = type)) + 
  theme_classic()+
  geom_line()



##########################################
# Run full model.
##########################################

user <- Sys.info()["user"]
store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"

jobname <-  "niv_sigma"

setwd(paste0(store_directory , "simsurvextrap_", jobname, "/"))
print(getwd())

scenarios <- read_csv("scenarios_niv_PFS.csv")
## View(scenarios)

scenarios <-  scenarios %>%
  filter(stan_fit_method == "mcmc")

for(i in scenarios$scenario_id){

  temp_sigma_draws <- readRDS(paste0("scen", i, "_sigma_draws.rds")) 
  
  if(i == 1){
    sigma_draws_all <- temp_sigma_draws
  } else{
    sigma_draws_all <- bind_rows(sigma_draws_all,
                                 temp_sigma_draws)
  }
  print(i)
  #sigma_draws 
}


saveRDS(sigma_draws_all, "sigma_draws_all.rds")


#####################################################
# Derive prior and posterior densities.
#####################################################

nsim <- 50

for(i in scenarios$scenario_id){
  #i <- 5
  for(j in 1:nsim){
  #j <- 1
  temp_dens <- density(sigma_draws_all %>%
                         filter(isim == j, scenario_id == i) %>%
                         pull(sigma), 
                       bw = 0.10)
  
  test_dens <- sigma_draws_all %>%
    filter(isim == j, scenario_id == i)
  
  #summary(test_dens)
  #  plot(temp_dens)
  
  temp_dens_df <- tibble(x = temp_dens$x, density = temp_dens$y) %>%
    mutate(scenario_id = i) %>%
    mutate(isim = j)
  
  
  if(j == 1){
    comb_dens_df <- temp_dens_df
  } else{
    comb_dens_df <- rbind(comb_dens_df, temp_dens_df)
  }
  
  
  }
  
  saveRDS(comb_dens_df, paste0("scen", i, "_sigma_draws_density.rds")) 
  
  print(i)
  #temp_dens_df
  #plot(temp_dens)

}


for(i in scenarios$scenario_id){
  
  #i <- 1
  
  post_df <- readRDS(paste0("scen", i, "_sigma_draws_density.rds")) 
  

  if(i == 1){
    comb_post_df <-  post_df
  } else{
    comb_post_df <- rbind(comb_post_df,  post_df)
  }
  
  
}
saveRDS(comb_post_df, paste0("sigma_post_density_all.rds")) 



for(i in scenarios$scenario_id){
  
  prior_df <- tibble(x = c(seq(from = 0, to = 0.5, length.out = 500),
                           seq(from = 0.5, to = 20, length.out = 500))
                     ) %>%
    mutate(density = dgamma(x = x, 
                         shape = scenarios$prior_hsd_shape[i], 
                         rate = scenarios$prior_hsd_rate[i])) %>%
    mutate(scenario_id = i) %>%
    mutate(isim = 0)
  
  saveRDS(prior_df, paste0("scen", i, "_sigma_prior_density.rds")) 
  
  print(i)
  #temp_dens_df
  #plot(temp_dens)
  
  if(i == 1){
    comb_prior_df <- prior_df
  } else{
    comb_prior_df <- rbind(comb_prior_df, prior_df)
  }
  
  
}

saveRDS(comb_prior_df, paste0("sigma_prior_density_all.rds")) 


#####################################################
# Plot prior posterior.
#####################################################


store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/"

jobname <-  "cetux_sigma"


setwd(paste0(store_directory , "simsurvextrap_", jobname, "/"))

if(jobname == "cetux_sigma"){
  scenarios <- read_csv("scenarios_cetux_OS.csv")
}else{
  scenarios <- read_csv("scenarios_niv_PFS.csv")
}


sigma_post_all <- readRDS("sigma_post_density_all.rds")
sigma_prior_all <- readRDS("sigma_prior_density_all.rds")

# head(sigma_post_all)
# head(sigma_prior_all)
# summary(sigma_post_all)
# summary(sigma_prior_all) 
# View(sigma_prior_all)

sigma_df <- rbind(sigma_post_all,
                  sigma_prior_all)

#summary(sigma_df)
#names(scenarios)
#summary(as.factor(scenarios$smooth_model))
#summary(scenarios$scenario_id)

for(smooth_model_type in c("random_walk", "exchangeable")){
  
  if(smooth_model_type == "random_walk") smooth_model_ext <- "randwalk"
  if(smooth_model_type == "exchangeable") smooth_model_ext <- "exch"
  
  scen_df <- scenarios %>%
    mutate(scenario_id = as.character(scenario_id)) %>%
    filter(mspline_df %in% c(3,6, 10),
           prior_hsd_rate %in% c(1,5,20),
           stan_fit_method %in% "mcmc",
           smooth_model %in% smooth_model_type,
           hscale_default %in% TRUE,
           mspline_bsmooth %in% TRUE) %>%
    mutate(label_df = paste0("`df=`*", mspline_df)) %>%
    mutate(label_hsd = paste0("sigma*","`~`*`G`*`amma`*`(`*", 
                              prior_hsd_shape, "*`,`*" , 
                              prior_hsd_rate, "*`)`")) %>% 
    mutate(label_scen = paste0("`df=`*", mspline_df, "*`,`*", "~sigma*","`~`*`G`*`amma`*`(`*", 
                               prior_hsd_shape, "*`,`*" , 
                               prior_hsd_rate, "*`)`" )) 
  
  plot_prior_posterior <- sigma_df %>%
    mutate(scenario_id = as.character(scenario_id)) %>%
    inner_join(scen_df, by = "scenario_id") %>%
    arrange(prior_hsd_rate) %>%
    arrange(mspline_df) %>%
    mutate(label_df = fct_inorder(label_df)) %>%
    mutate(label_hsd = fct_inorder(label_hsd)) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_colour = if_else(isim  == 0, 1, 0),
           line_colour = as.factor(line_colour)) %>%
    ggplot(aes(x = x, y= density, alpha = line_alpha, colour = line_colour, 
               group = isim, linewidth = line_width))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 6),
          legend.box.margin =  margin(0, 0, 0,0),
          legend.box.spacing = unit(0.2, "pt"),
          legend.text = element_markdown(size = 10)) +
     geom_line()+
    scale_linewidth(NULL, range = c(0.4,1))+
    scale_alpha(range = c(0.1,1))+
    scale_color_manual(NULL, values = c("#830051", "gray30"),
                       labels = c(
                         "<b>&sigma; posterior</b><br><span style='font-weight:400;'>(50 simulations)</span>",
                         "<b>&sigma; prior</b>"
                       ))+
    facet_wrap(~label_df+label_hsd, ncol = 3, labeller = label_parsed)+
    scale_x_continuous(NULL, limit = c(0,5.2), breaks = 1:5)+
    scale_y_continuous("Density", limit = c(0,4.3), breaks = 1:4)+
    guides(linewidth = "none", alpha = "none")  +
    guides(colour = guide_legend(override.aes = list(linewidth = c(0.6,1))))
 # plot_prior_posterior
 # assign(paste0("sigma_", smooth_model_ext), plot_prior_posterior)
  saveRDS(plot_prior_posterior, paste0("sigma_", smooth_model_ext, "_prior_posterior_figure.rds"))
}

store_dir_cetux <- paste0(store_directory , "simsurvextrap_cetux_sigma/")
store_dir_niv <- paste0(store_directory , "simsurvextrap_niv_sigma/")

#setwd()
sigma_cetux_rw <- readRDS(paste0(store_dir_cetux, "sigma_randwalk_prior_posterior_figure.rds"))
sigma_niv_rw <- readRDS(paste0(store_dir_niv, "sigma_randwalk_prior_posterior_figure.rds"))
sigma_cetux_exch <- readRDS(paste0(store_dir_cetux, "sigma_exch_prior_posterior_figure.rds"))
sigma_niv_exch <- readRDS(paste0(store_dir_niv, "sigma_exch_prior_posterior_figure.rds"))


figure_all_rw <- plot_grid(
  sigma_cetux_rw  +
    theme(
      plot.title = element_markdown(size = 10, face = "bold"),
      plot.title.position = "plot"
    ) +
    labs(title = "(a) Cetuximab OS, <b>&sigma;</b> prior and posterior"),
  sigma_niv_rw+ 
    theme(
      plot.title = element_markdown(size = 10, face = "bold"),
      plot.title.position = "plot"
    ) +
    labs(title = "(b) Nivolumab PFS, <b>&sigma;</b> prior and posterior"), 
  ncol = 1,
  align = "v")

figure_all_exch <- plot_grid(
  sigma_cetux_exch +
    theme(
      plot.title = element_markdown(size = 10, face = "bold"),
      plot.title.position = "plot"
    ) +
    labs(title = "(a) Cetuximab OS, <b>&sigma;</b> prior and posterior"),
  
  sigma_niv_exch +
    theme(
      plot.title = element_markdown(size = 10, face = "bold"),
      plot.title.position = "plot"
    ) +
    labs(title = "(b) Nivolumab PFS, <b>&sigma;</b> prior and posterior"),
  
  ncol = 1,
  align = "v"
)


tiff(file = paste0(store_directory, "figures/figure_sigma_priors_randwalk.tiff"),   
     width = 6.5, 
     height = 8,
     units = 'in',  
     res = 1200, 
     compression = "lzw")
figure_all_rw
dev.off()

tiff(file = paste0(store_directory, "figures/figure_sigma_priors_exch.tiff"),   
     width = 6.5, 
     height = 8,
     units = 'in',  
     res = 1200, 
     compression = "lzw")
figure_all_exch
dev.off()

