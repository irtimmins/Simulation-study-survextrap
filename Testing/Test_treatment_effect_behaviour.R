#####


library(pracma)
library(survival)
library(flexsurv)
library(dplyr)
library(tidyr)
library(ggplot2)

#   # 0.8*(exp(-(t-4)))*(-0.02-0.02*t+0.029*t^2+0.018*t^3-0.0059*t^4)*trt

cetux <- readRDS("Data/cetuximab_OS.rds") 
surv_df <- cetux[cetux$treat=="Cetuximab",]
# end of trial follow-up time for rmst computation
maxT_data <- max(surv_df$years)
# specify number of knots for spline model
k_true <-  3
# derive true model using flexsurv
true_mod <- flexsurv::flexsurvspline(Surv(years, d) ~ 1, data=surv_df, k=k_true)
plot(true_mod, type = "hazard")
t_val <- 2
GL <- gaussLegendre(n = 100, 0, t_val)
GL_times <- GL$x
GL_weights <- GL$w
beta2 <- c(log(0.5)/(1-tanh(-1.2)), 0.8, -1.2)
beta_true <- beta2

gamma <- true_mod$coefficients
knots <- true_mod$knots 
beta <-  beta_true 
#trt <- 0

a <- log(0.5); b <- log(1.8)
ks <- 5; kd <- 10

s <- function(t) 0.5*(1 + tanh(ks*(t - 1)))
d <- function(t) 1 - 1/(1 + exp(-kd*(t - 2.5)))

beta4 <- c(5, 10, log(0.5), log(1.8))

haz <- function(t, gamma, knots, beta, trt) {
 
  s <- function(t) 0.5*(1 + tanh(beta[1]*(t - 1)))
  d <- function(t) 1 - 1/(1 + exp(-beta[2]*(t - 2.5)))
  
  loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+ 
  (beta[3] + (beta[4] - beta[3]) * s(t)) * d(t)*trt
  
  return(exp(loghaz))
} 


log(hsurvspline(x = 0.1, gamma = gamma, knots = knots))
exp(log(hsurvspline(x = 0.1, gamma = gamma, knots = knots)))
#t <- 0
#log(0.2)*exp(-0.05*t)+4.5*(exp(-(t-2)^2/(2*log(2)^2))-exp(-2/(log(2)^2))*exp(-0.05*t))
haz(t = 0.1, gamma = gamma, knots = knots, beta = beta_true, trt = 1)/
  haz(t = 0.1, gamma = gamma, knots = knots, beta = beta_true, trt = 0)

haz(t = 0.1, gamma = gamma, knots = knots, beta = beta_true, trt = 0)

haz_control <- function(t){
  haz_out <- haz(t, gamma = gamma, knots = knots, beta = beta4, trt = 0)
  return(haz_out)
}
haz_control(2.5)

haz_active <- function(t){
  haz_out <- haz(t, gamma = gamma, knots = knots, beta = beta4, trt = 1)
  return(haz_out)
}
haz_control(10)
haz_active(c(10))
haz_control(0.1)/haz_active(0.1)

haz_control(0.05)

survival_control <- function(t){
  
  GL <- gaussLegendre(n = 100, 0, t)
  GL_times <- GL$x
  GL_weights <- GL$w
  survival_out <- exp(-sum(haz_control(GL_times)*GL_weights))
  
  return(survival_out)
  
}

survival_active <- function(t){
  
  GL <- gaussLegendre(n = 100, 0, t)
  GL_times <- GL$x
  GL_weights <- GL$w
  survival_out <- exp(-sum(haz_active(GL_times)*GL_weights))
  
  return(survival_out)
  
}


haz_test <- tibble(t = seq(from = 0, to = 5, length.out = 100)) %>%
  mutate(control = 0, active = 0) %>%
  mutate(control = haz_control(t),
         active = haz_active(t))

haz_test %>%
  mutate(hazard_ratio = active/control) %>%
  #pivot_longer(cols = c(control, active), values_to = "hazard", names_to = "arm") %>%
  ggplot(aes(x = t, y= hazard_ratio))+
  theme_classic()+
  geom_line()+
  scale_y_continuous(limits = c(0,2.5))+
  geom_hline(yintercept = 1)


data_test <- tibble(t = seq(from = 0, to = 5, length.out = 100)) %>%
  mutate(control = 0, active = 0)

for(i in 1:100){
data_test$control[i] <- survival_control(data_test$t[i])
data_test$active[i] <- survival_active(data_test$t[i])
}

data_test %>%
  pivot_longer(cols = c(control, active), values_to = "survival", names_to = "arm") %>%
  ggplot(aes(x = t, y= survival, colour = arm ))+
  theme_classic()+
  geom_line()+
  scale_y_continuous(limits = c(0,1))

#survival_control(2)



