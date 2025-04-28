## Simulate from the "true" survival model

#' @param model_class class of true model e.g. flexsurvreg
#' @param true_mod true model 
#' @param N Number of patients in trial
#' @param nsim number of simulated datasets
#' @param maxT maximum event time
#' @param seed Random number seed


### for review

sim_dgm <- function(model_class = "flexsurv", true_mod, N, nsim, maxT, seed = NULL){
  
  set.seed(seed)
  
  
  # Create a new simulated baseline characteristics dataset of size N
  
  # may need to adapt code for different classes 
  # flexsurvspline, flexsurvreg, flexsurvmix, survextrap...
  # create individual helper functions for each.
  
  x <- data.frame(id=1:N)
  
  # Could also add a censoring distribution. If this is done will need to change where seed is specified.
  x$censtime <- maxT
  
  # Simulate events
  sim_data <- simulate(true_mod, nsim = nsim, newdata = x,
                       seed = seed,
                       tidy=T, censtime = x$censtime)
  
  
  return(sim_data=sim_data)
  
  
}



sim_dgm_mix <- function(mix_params = NULL, gpm_params = NULL, age_dist = NULL, 
                        N, nsim, maxT, seed = NULL){
  
  set.seed(seed)
  
  gpm_include_true <- !(anyNA(gpm_params) | is.null(gpm_params))
  
  if(gpm_include_true) mix_params <- c(mix_params, gpm_params)
  
  
  x_mix <- tibble(id = 1:(N*nsim)) %>%
    mutate(z = rnorm(n = n(), mean = 0, sd = 1),
           age = rnorm(N*nsim, mean = age_dist$mean, sd = age_dist$sd))
  
  cumhaz.mix <- function(t, x, betas){
    
    cumh <-  -exp(x[["z"]]*betas[["frailty"]])*
      log(
        betas[["pmix"]]*
          exp(-betas[["lambda1"]]*t^betas[["gamma1"]])+
          (1-betas[["pmix"]])*
          exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
      )
    
    res <- cumh
    
    res
  }
  
  
  cumhaz.mix.gpm <- function(t, x, betas){
    
    cumh <-  -exp(x[["z"]]*betas[["frailty"]])*
      log(
        betas[["pmix"]]*
          exp(-betas[["lambda1"]]*t^betas[["gamma1"]])+
          (1-betas[["pmix"]])*
          exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
      )
    
    cumh.gpm <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*(t+x[["age"]]))-1)
    
    cumh.gpm.cens <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*x[["age"]])-1)
    
    res <- cumh + cumh.gpm - cumh.gpm.cens
    
    res
    
  }
  
  
  sim_data <- simsurv(betas = mix_params,
                      x = x_mix,
                      interval = c(10^(-20), 10^20), 
                      cumhazard = ifelse(gpm_include_true, cumhaz.mix.gpm, cumhaz.mix),
                      maxt = maxT)

  sim_data <- sim_data %>%
    select(-id) %>%
    mutate(censtime = maxT) %>%
    rename(time = eventtime,
           event = status) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim)
      ) %>%
    relocate(id, censtime, i, time, event)
  
  return(sim_data=sim_data)
}


# Function to simulate mixture-Weibull for disease-specific mortality and a Gompertz for other-cause mortality
# Combines to get all-cause mortality
#' @param mix_params parameters of mixture Weibull model for disease-specific mortality
#' @param gpm_params parameters of GPM mortality
#' @param age_dist mean and SD of ages
#' @param N Number of patients in trial
#' @param nsim number of simulated datasets
#' @param maxT maximum event time
#' @param seed Random number seed
#' @param trt_beta Maximum treatment effect during stable treatment efficacy period t2<t<t3
#' @param trt_timepoints Time points that define 1) period of delayed treatment effect 0<t<t1, 
#' 2) treatment efficacy increasing t1<t<t2, 3) treatment effect constant t2<t<t3, 4) treatment effect waning t3<t<t4,
#' 5) treatment effect null t4<t
sim_dgm_mix_haz <- function(mix_params = NULL, gpm_params = NULL, age_dist = NULL, 
                              N, nsim, maxT, seed = NULL, 
                              trt_beta = NULL, trt_timepoints = NULL){
  
  # Set seed
  set.seed(seed)
  
  gpm_include_true <- !(anyNA(gpm_params) | is.null(gpm_params))
  
    if(gpm_include_true) mix_params <- c(mix_params, gpm_params)
  
  # generate covariate distributions for N*nsim patients
  # N is number of patients in a trial, so N/2 per arm
  x_mix <- tibble(id = 1:(N*nsim)) %>%
    mutate(z = rnorm(n = n(), mean = 0, sd = 1),
           age = rnorm(n = n(), mean = age_dist$mean, sd = age_dist$sd)) %>%
    mutate(trt = rep_len(0:1, n()))
  
  
  haz_trt <- function(t){
    
    if(t <= trt_timepoints[1]) return(0)
    if(t > trt_timepoints[1] & t <= trt_timepoints[2]) return(trt_beta*(t-trt_timepoints[1])/(trt_timepoints[2]-trt_timepoints[1]))
    if(t > trt_timepoints[2] & t <= trt_timepoints[3]) return(trt_beta) 
    if(t > trt_timepoints[3] & t <= trt_timepoints[4]) return(trt_beta*(1-(t-trt_timepoints[3])/(trt_timepoints[4]-trt_timepoints[3])))
    if(t > trt_timepoints[4]) return(0)
    
  }
  
  haz.mix <- function(t, x, betas){
    
    pdf1 <- betas[["lambda1"]]*betas[["gamma1"]]*(t^(betas[["gamma1"]]-1))*exp(-betas[["lambda1"]]*t^betas[["gamma1"]])
    pdf2 <- betas[["lambda2"]]*betas[["gamma2"]]*(t^(betas[["gamma2"]]-1))*exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
    S1 <- exp(-betas[["lambda1"]]*t^betas[["gamma1"]])
    S2 <- exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
    
    basehaz <- (betas[["pmix"]]*pdf1 + (1-betas[["pmix"]])*pdf2) /
      (betas[["pmix"]]*S1+ (1-betas[["pmix"]])*S2)
    
    haz <- basehaz * exp(x[["trt"]]*haz_trt(t))
    
    res <- haz
    
    res
  }
  
  
  cumhaz.gpm <- function(t, x, betas){

    cumh.gpm <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*(t+x[["age"]]))-1)
    
    cumh.gpm.cens <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*x[["age"]])-1)
    
    res <- cumh.gpm - cumh.gpm.cens
    
    res
    
  }
  
  
  sim_data.cancer <- simsurv(betas = mix_params,
                      x = x_mix,
                      interval = c(10^(-20), maxT+1), # speed up considerably by limiting interval to maxT+1
                      hazard = haz.mix,
                      maxt = maxT)
  
  sim_data.gpm <- simsurv(betas = mix_params,
                      x = x_mix,
                      interval = c(10^(-20), maxT+1), # speed up considerably by limiting interval to maxT+1
                      cumhazard = cumhaz.gpm,
                      maxt = maxT)
  
  #nrow(sim_data)
  sim_data <- sim_data.cancer %>%
    select(-id) %>%
    mutate(censtime = maxT) %>%
    rename(time_cancer = eventtime,
           event_cancer = status) %>%
    mutate(time_gpm = sim_data.gpm[["eventtime"]],
           event_gpm = sim_data.gpm[["status"]]) %>%
    mutate(time = pmin(time_cancer, time_gpm),
           event = ifelse(time==time_cancer, event_cancer, event_gpm)) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim),
      x_mix["trt"],
      .before = 1
    ) 
    
  return(sim_data=sim_data)
}




###################################################
# test functions
###################################################


sim_dgm_mix_try <- function(mix_params = NULL, gpm_params = NULL, age_dist = NULL, 
                        N, nsim, maxT, seed = NULL){
  
  set.seed(seed)
  
  gpm_include_true <- !(anyNA(gpm_params) | is.null(gpm_params))
  
  if(gpm_include_true) mix_params <- c(mix_params, gpm_params)
  

  cumhaz.mix <- function(t, x, betas){
    
    cumh <-  -exp(x[["z"]]*betas[["frailty"]])*
      log(
        betas[["pmix"]]*
          exp(-betas[["lambda1"]]*t^betas[["gamma1"]])+
          (1-betas[["pmix"]])*
          exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
      )
    
    res <- cumh
    
    res
  }
  
  
  cumhaz.mix.gpm <- function(t, x, betas){
    
    cumh <-  -exp(x[["z"]]*betas[["frailty"]])*
      log(
        betas[["pmix"]]*
          exp(-betas[["lambda1"]]*t^betas[["gamma1"]])+
          (1-betas[["pmix"]])*
          exp(-betas[["lambda2"]]*t^betas[["gamma2"]])
      )
    
    cumh.gpm <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*(t+x[["age"]]))-1)
    
    cumh.gpm.cens <- betas[["lambda_gompertz"]]*(betas[["gamma_gompertz"]]^(-1))*(exp(betas[["gamma_gompertz"]]*x[["age"]])-1)
    
    res <- cumh + cumh.gpm - cumh.gpm.cens
    
    res
    
  }
  
  
  sim_data <- NULL
  
  # simulate data iteratively in K blocks of up to 10,000.
  
  Kq <- (N*nsim) %/% 10000
  Kr <- (N*nsim) %% 10000
  
  # print(paste0("Kq ", Kq))
  # print(paste0("Kr ", Kr))
  
  K <- if_else(Kr == 0, Kq,  Kq + 1 )
  
  #print(K)
  
  for(i in 1:K){
    
    sim_success <- FALSE
    count_tries <- 0
    
    nblock <- if_else(i < K | Kr == 0, 10000,  Kr)
    
    # print(paste0("nblock ", nblock))
    
    x_mix <- tibble(id = 1:nblock) %>%
      mutate(z = rnorm(n = nblock, mean = 0, sd = 1),
             age = rnorm(nblock, mean = age_dist$mean, sd = age_dist$sd))
   
    #print(paste0("x_mix length ", nrow(x_mix)))
    
    while(sim_success == FALSE & count_tries <= 50){  
      temp_data <- try(simsurv(betas = mix_params,
                                           x = x_mix,
                                           interval = c(10^(-20), 10^20), 
                                           cumhazard = ifelse(gpm_include_true, cumhaz.mix.gpm, cumhaz.mix),
                                           maxt = maxT))
      
      sim_success <- (class(temp_data) != "try-error")
      count_tries <- 1 + count_tries
      #print(paste0("count ",count_tries))
    }
    
    if(i == 1) {
      sim_data <- temp_data
    } else {
      sim_data <- rbind(sim_data, temp_data)
    }
    
    #print(paste0(round(100*nrow(sim_data)/N), "%"))
    
  }
  
  #print(count_tries)
  #print(nrow(sim_data))
  
  sim_data <- sim_data %>%
    select(-id) %>%
    mutate(censtime = maxT) %>%
    rename(time = eventtime,
           event = status) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim)
    ) %>%
    relocate(id, censtime, i, time, event)
  
  return(sim_data=sim_data)
}
















sim_dgm_mix2 <- function(mix_params = NULL, gpm_params = NULL, age_dist = NULL, 
                        N, nsim, maxT, seed = NULL){
  
  
  set.seed(seed)
  
  
  x_mix <- tibble(id = 1:(N*nsim)) %>%
    mutate(z = rnorm(n = n(), mean = 0, sd = 1),
           age = rnorm(N*nsim, mean = age_dist$mean, sd = age_dist$sd))
  
  
  gpm_include_true <- !(anyNA(gpm_params) | is.null(gpm_params))
  
  if(gpm_include_true) mix_params <- c(mix_params, gpm_params)
  

  sim_data <- simsurv(lambdas = mix_params[c("lambda1", "lambda2")],
                      gammas = mix_params[c("gamma1", "gamma2")],
                      mixture = TRUE,
                      pmix = mix_params["pmix"],
                      x = x_mix,
                      interval = c(10^(-20), 10^20), 
                      maxt = maxT)
  
  sim_data <- sim_data %>%
    select(-id) %>%
    mutate(censtime = maxT) %>%
    rename(time = eventtime,
           event = status) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim)
    ) %>%
    relocate(id, censtime, i, time, event)
  
  return(sim_data=sim_data)
}


sim_dgm_mix3 <- function(mix_params = NULL, N, nsim, maxT, seed = NULL){
  
  
  set.seed(seed)

  x_mix <- tibble(id = 1:(N*nsim)) %>%
    mutate(z = rnorm(n = n(), mean = 0, sd = 1))
  
  
  sim_data1 <- simsurv(lambdas = mix_params["lambda1"],
                      gammas = mix_params["gamma1"],
                      dist = "weibull",
                      mixture = FALSE,
                      x = x_mix,
                      interval = c(10^(-20), 10^20), 
                      maxt = maxT)
  
  sim_data2 <- simsurv(lambdas = mix_params["lambda2"],
                       gammas = mix_params["gamma2"],
                       dist = "weibull",
                       mixture = FALSE,
                       x = x_mix,
                       interval = c(10^(-20), 10^20), 
                       maxt = maxT)
  
  sim_data <- sim_data1 %>%
    select(-id) %>%
    mutate(censtime = maxT) %>%
    rename(time1 = eventtime,
           event1 = status) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim)
    ) %>%
    relocate(id, censtime, i, time1, event1)  %>%
    add_column(time2 = sim_data2[["eventtime"]],
               event2 = sim_data2[["status"]] ) %>%
    mutate(gen_unif = runif(n = n()),
           gen_unif = 1*(gen_unif <= mix_params["pmix"])) %>%
    mutate(time = 0,
           event = 0) %>%
    mutate(time = case_when(gen_unif == 1 ~ time1,
                            gen_unif == 0 ~ time2),
           event = case_when(gen_unif == 1 ~ event1,
                            gen_unif == 0 ~ event2)) %>%
    select(id, censtime, i, time, event)
  

  return(sim_data=sim_data)
}




sim_dgm_mix3_try <- function(mix_params = NULL, N, nsim, maxT, seed = NULL){
  
  
  set.seed(seed)
  
  
  
  sim_data <- NULL
  
  # simulate data iteratively in K blocks of up to 10,000.
  
  Kq <- (N*nsim) %/% 10000
  Kr <- (N*nsim) %% 10000
  
  #print(paste0("Kq ", Kq))
  #print(paste0("Kr ", Kr))
  
  K <- if_else(Kr == 0, Kq,  Kq + 1 )
  
  #print(K)
  
  for(i in 1:K){
    
    sim_success <- FALSE
    count_tries <- 0
    
    nblock <- if_else(i < K | Kr == 0, 10000,  Kr)
    
    #print(paste0("nblock ", nblock))
    
    x_mix <- tibble(id = 1:nblock) 
    
    #print(paste0("x_mix length ", nrow(x_mix)))
    
    while(sim_success == FALSE & count_tries <= 50){  
      temp_data1 <- try(simsurv(lambdas = mix_params["lambda1"],
                               gammas = mix_params["gamma1"],
                               dist = "weibull",
                               mixture = FALSE,
                               x = x_mix,
                               interval = c(10^(-20), 10^20), 
                               maxt = maxT))
      
      temp_data2 <- try(simsurv(lambdas = mix_params["lambda2"],
                                gammas = mix_params["gamma2"],
                                dist = "weibull",
                                mixture = FALSE,
                                x = x_mix,
                                interval = c(10^(-20), 10^20), 
                                maxt = maxT))
      
      sim_success <- (class(temp_data1) != "try-error") & (class(temp_data2) != "try-error")
      count_tries <- 1 + count_tries
     # print(paste0("count ",count_tries))
    }
    

temp_data <- temp_data1 %>%
    mutate(censtime = maxT) %>%
    rename(time1 = eventtime,
           event1 = status) %>%
    add_column(time2 = temp_data2[["eventtime"]],
               event2 = temp_data2[["status"]] ) %>%
    mutate(gen_unif = runif(n = n()),
           gen_unif = 1*(gen_unif <= mix_params["pmix"])) %>%
    mutate(time = 0,
           event = 0) %>%
    mutate(time = case_when(gen_unif == 1 ~ time1,
                            gen_unif == 0 ~ time2),
           event = case_when(gen_unif == 1 ~ event1,
                             gen_unif == 0 ~ event2)) 
    
    
    if(i == 1) {
      sim_data <- temp_data
    } else {
      sim_data <- rbind(sim_data, temp_data)
    }
    
    #print(paste0(round(100*nrow(sim_data)/N), "%"))
    
  }
  
  sim_data <- sim_data %>%
    select(-id) %>%
    add_column(
      expand_grid(id = 1:N, i = 1:nsim)
    ) %>%
    select(id, censtime, i, time, event)  
  
  
  return(sim_data=sim_data)
}

