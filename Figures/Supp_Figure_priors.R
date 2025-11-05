library(survextrap)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(cowplot)
library(purrr)

store_directory <- "C:/Users/LSHIT9/OneDrive - London School of Hygiene and Tropical Medicine/Documents/Survival extrapolation/Paper 1/Data/figures"

surv_df <- readRDS("Data/cetuximab_OS_simulated_example.rds")

for(k in 1:3){
  
  
  prior_choice <- c(1, 5, 20)[k]
  
 # mspline_evenly_spaced = list(knots = seq(from = 0.5, to = 5, length.out = 9))
  
  nd_mod <- survextrap(Surv(time, event) ~ 1, data=surv_df, smooth_model =  "random_walk", fit_method = "opt",
                       prior_hsd = p_gamma(2, prior_choice)) #, mspline = mspline_evenly_spaced)
  #nd_mod$mspline$knots
  
  
  #plot(nd_mod)
  #print_priors(nd_mod)
  
  
  
  ###############################################
  # Hazard trajectories.
  ###############################################
  
  set.seed(15151)
  haz_sim <- nd_mod$prior_sample$haz(nsim=200)
  plot_a <- ggplot(haz_sim, aes(x=time, y=haz, group=rep)) + 
    theme_classic()+
    geom_line( colour = "gray50", alpha = 0.7) +
    scale_y_continuous("Hazard", limits=c(0,54))+
    scale_x_continuous("Time (years)", breaks = 0:5, limits = c(0,5))
  
  plot_a
  
  ###############################################
  # p_i distributions.
  ###############################################
  
  
  #nd_mod$stanfit$par
  set.seed(1251)
  mod_samples <- nd_mod$prior_sample$sample(nsim = 100000)
  #mod_samples
  mod_samples_coef <- mod_samples$coefs
  #names(mod_samples)
  #mod_samples$coefs
  
  #sum(mod_samples_coef[2,1:10])
  
  plot_b <- mod_samples_coef %>%
    as_tibble() %>%
    mutate(iter = row_number()) %>%
    rename_with(
      .fn = ~ paste0("italic(p)[", seq_along(.x), "]"),   # create new names p1, p2, ...
      .cols = starts_with("v")              # select columns that start with 'v'
    ) %>%
    pivot_longer(names_to = "weights", cols = starts_with("italic"))%>%
    mutate(
      id_num = readr::parse_number(weights),
      weights = factor(weights, levels = paste0("italic(p)[", sort(unique(id_num)), "]"))) %>%
    ggplot()+
    theme_classic()+
    geom_density(aes(x = value), colour = "gray30", alpha = 0.7)+
    facet_wrap(~weights, labeller = label_parsed)+
    ylab("Density")+
    ylim(c(0,100))+
    scale_x_continuous(NULL, breaks = c(0, 0.5, 1), limit = c(0,1))
  
  plot_b
  
  ###############################################
  # p1*b1(t) trajectories.
  ###############################################
  
  const_haz_prior <- mspline_constant_coefs(
    nd_mod$mspline
  )
  const_haz_prior
  
  p_vec_samples <- mod_samples_coef %>%
    as_tibble() %>%
    slice(1:100)
  
  
  for(i in 1:25){
    
    temp_p_vec <- p_vec_samples[i,] %>% 
      as.numeric()
    
    for(j in 1:10){
      
      scale_factor <- temp_p_vec[j]
      
      temp_data <- 
        mspline_plotdata(knots = nd_mod$mspline$knots,
                         coefs =  diag(nd_mod$mspline$df)[j,],
                         scale = scale_factor) %>%
        as_tibble() %>%
        filter(term == j) %>%
        mutate(
          id_num = term,
          term_label = paste0("italic(p)[", id_num,
                              "]", "*italic(b)[", id_num, "](t)")) %>%
        mutate(term = factor(term, levels = term, labels = term_label)) %>%
        select(time, term, haz, id_num) %>%
        mutate(iter = i)
      
      if(j == 1) prior_data_i <- temp_data
      if(j != 1) prior_data_i <- bind_rows( prior_data_i,
                                            temp_data)
    }
    print(i) 
    
    if(i == 1) prior_data_all <- prior_data_i
    if(i != 1) prior_data_all <- bind_rows(prior_data_all,
                                           prior_data_i)
    
  }
  
  plot_c <- prior_data_all %>%
    ggplot(aes(x = time, y = haz, group = iter))+
    theme_classic()+
    geom_line(colour = "gray50", alpha = 0.7)+
    facet_wrap(~term, labeller = label_parsed)+
    scale_y_continuous(NULL, limits = c(0, 1.05), breaks = c(0,1))+
    scale_x_continuous("Time (years)", limits = c(0,5))
  
  plot_label <- c("a", "b", "c")[i]
  #plot_c
  plot_all <- plot_grid(plot_a+
                          theme(plot.title = element_text(vjust =0.2,
                                                          size=10, 
                                                          face="bold", lineheight = 0.2),
                                plot.title.position = "plot")+
                          labs(
                            title = bquote(
                              bold(
                                atop(
                                  "Priors and hazard trajectories for " ~ 
                                    sigma ~ "~ Gamma(2," * .(as.character(prior_choice))* ")        ",
                                  "(a) Samples from the prior distribution for the hazard    "
                                )
                              )
                            )
                          ),
                        plot_b+
                          theme(plot.title = element_text(vjust =0.2,
                                                          size=10, 
                                                          face="bold"),
                                plot.title.position = "plot")+
                          labs(title = expression(bold("(b) Prior distribution of the M-spline weights") ~ bolditalic(p)[i])),
                        plot_c+
                          theme(plot.title = element_text(vjust =0.2,
                                                          size=10, 
                                                          face="bold"),
                                plot.title.position = "plot")+
                          labs(title = expression(bold("(c) Samples from the prior distribution of") ~ bolditalic(p)[i] * bolditalic(b)[i] *
                                                    "(" * bold(t) * ")")),
                        ncol = 1,
                        rel_heights = c(0.28, 0.36, 0.36))
  
  
  plot_all
  #i<- 1
  
  
  pdf(file = paste0(store_directory, "/figure_prior_", k, ".pdf"),   
      width = 5.4, 
      height = 8.4) 
  print(plot_all) 
  dev.off()
  
  
  tiff(file = paste0(store_directory, "/figure_prior_", k, ".tiff"),   
       width = 5.4, 
       height = 8.4,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all) 
  dev.off()
  
  
  #?plot_grid()
  #plot_all
  
}

