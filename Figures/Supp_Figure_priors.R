library(survextrap)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(purrr)


surv_df <- readRDS("Data/cetuximab_OS.rds")
nd_mod <- survextrap(Surv(years, d) ~ 1, data=surv_df, fit_method="opt", iter = 2000)
#print_priors(nd_mod)

#nd_mod$stanfit$par

mod_samples <- nd_mod$prior_sample$sample(nsim = 1000)
mod_samples
mod_samples_coef <- mod_samples$coefs
#names(mod_samples)
#mod_samples$coefs

#sum(mod_samples_coef[2,1:10])

plot1 <- mod_samples_coef %>%
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
    geom_density(aes(x = value))+
    facet_wrap(~weights, labeller = label_parsed)+
  xlab(NULL)+
  ylab("Density")+
  xlim(c(0,1))+
  ylim(c(0,20))+
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))


######################
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

prior_data_all %>%
ggplot(aes(x = time, y = haz, group = iter))+
  theme_classic()+
  geom_line()+
  facet_wrap(~term, labeller = label_parsed)+
  scale_y_continuous(NULL)+
  scale_x_continuous("Time (years)", limits = c(0,5))
  






