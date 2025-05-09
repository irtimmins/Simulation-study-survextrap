## Evaluate estimands for true/fitted models for single arm models.

#' @param sim_data_chr Name of dataset to extract simulated clinical trial data. Dataset includes N clinical trials.
#' @param isim Iteration to extract from simulated clinical trial data from isim in 1, ..., N.
#' @param prior_chr Character specifying object with prior list object for survextrap.
#' @param mspline_chr Character specifying object with m-spline list object for survextrap.
#' @param stan_fit_chr Character specifying Stan model fit method e.g. "mcmc" for survextrap.
#' @param chains_chr Character specifying object with number of MCMC chains for survextrap.
#' @param iter_chr Character specifying object with number of MCMC iterations for survextrap.
#' @param smooth_model_chr Character specifying object with smooth_model specification for survextrap.
#' @param save_file File location to store estimates.
#' @param seed Set seed for survextrap model fit.
#' @param num Redundant input used for post slurm checking.
#

fit_est_slurm <- function(sim_data_chr,
                          isim,
                          prior_chr,
                          mspline_chr,
                          stan_fit_chr,
                          chains_chr,
                          iter_chr,
                          smooth_model_chr,
                          save_file,
                          seed,
                          num){
  
  start_all <- Sys.time()
  
  set.seed(seed)
  
  sim_data <- get(sim_data_chr)
  prior <- get(prior_chr)
  mspline <- get(mspline_chr)
  stan_fit <- get(stan_fit_chr)
  smooth_model <- get(smooth_model_chr)
  
  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  start_fit <- Sys.time()
  
  
  if(stan_fit == "mcmc"){  
    chains <- get(chains_chr)
    iter <- get(iter_chr)
  } else {
    chains <- iter <- NA
  }
  
  start_fit <- Sys.time()
  
  fitted_mod <- survextrap(Surv(time, event) ~ 1, 
                           mspline = mspline, 
                           prior_hscale = prior$prior_hscale,
                           prior_hsd = prior$prior_hsd,
                           data = data,
                           fit_method=stan_fit$fit_method,
                           chains = chains,
                           iter = iter,
                           smooth_model = smooth_model) 
  
  end_fit <- Sys.time()
  
  diags_df <- diags(fitted_mod) %>% 
    rename(estimand=name) 
  
  start_est <- Sys.time()
  
  est_df <- mod_est(fitted_mod, estimand = NULL, t = t_vec, rmst_t = rmst_vec, nonprop = FALSE,
                    return_model = TRUE, save_file = NULL)
  
  end_est <- Sys.time()
  
  end_all <- Sys.time()
  
  
  timings <- tibble("estimand" = c("run_time", "fit_time", "est_time"),
                    "value" = c(as.numeric(difftime(end_all, start_all, units='mins')),
                                as.numeric(difftime(end_fit, start_fit, units='mins')),
                                as.numeric(difftime(end_est, start_est, units='mins'))
                    )  )
  
  
  res <- est_df %>% 
    full_join(diags_df, by=c("estimand","value")) %>%
    full_join(timings, by=c("estimand","value") )
  
  if (!is.null(save_file)) saveRDS(res, file=save_file)
  
}


## Evaluate estimands for true/fitted models for two arm trials.

#' @param sim_data_chr Name of dataset to extract simulated clinical trial data. Dataset includes N clinical trials.
#' @param isim Iteration to extract from simulated clinical trial data from isim in 1, ..., N.
#' @param prior_chr Character specifying object with prior list object for survextrap.
#' @param mspline_chr Character specifying object with m-spline list object for survextrap.
#' @param stan_fit_chr Character specifying Stan model fit method e.g. "mcmc" for survextrap.
#' @param chains_chr Character specifying object with number of MCMC chains for survextrap.
#' @param iter_chr Character specifying object with number of MCMC iterations for survextrap.
#' @param smooth_model_chr Character specifying object with smooth_model specification for survextrap.
#' @param save_file File location to store estimates.
#' @param seed Set seed for survextrap model fit.
#' @param num Redundant input used for post slurm checking.
#

fit_est_slurm_trt <- function(sim_data_chr,
                          isim,
                          nonprop,
                          prior_chr,
                          mspline_chr,
                          stan_fit_chr,
                          chains_chr,
                          iter_chr,
                          smooth_model_chr,
                          save_file,
                          seed,
                          num){
  
  start_all <- Sys.time()
  
  set.seed(seed)
  
  sim_data <- get(sim_data_chr)
  prior <- get(prior_chr)
  mspline <- get(mspline_chr)
  stan_fit <- get(stan_fit_chr)
  smooth_model <- get(smooth_model_chr)
  print(smooth_model)
  
  if(stan_fit == "mcmc"){  
    chains <- get(chains_chr)
    iter <- get(iter_chr)
  } else {
    chains <- iter <- NA
  }
  
  data <- sim_data %>%
    filter(i == isim) %>%
    select(-i)
  
  print(nrow(data))
  
  start_fit <- Sys.time()
  
  fitted_mod <- survextrap(Surv(time, event) ~ trt,
                           nonprop = nonprop,
                           mspline = mspline, 
                           prior_hscale = prior$prior_hscale,
                           prior_hsd = prior$prior_hsd,
                           prior_hrsd = prior$prior_hrsd,
                           data = data,
                           fit_method=stan_fit$fit_method,
                           chains = chains,
                           iter = iter,
                           smooth_model = smooth_model)  
  
  
  end_fit <- Sys.time()
  
  print(fitted_mod)
  
  diags_df <- diags(fitted_mod) %>% 
    rename(estimand=name) 
  
  start_est <- Sys.time()
  
  est_df <- mod_est_trt(fitted_mod, estimand = NULL, t = t_vec, rmst_t = rmst_vec, nonprop = FALSE,
                        return_model = TRUE, save_file = save_file)
  
  end_est <- Sys.time()
  
  end_all <- Sys.time()
  
  timings <- tibble("estimand" = c("run_time", "fit_time", "est_time"),
                    "value" = c(as.numeric(difftime(end_all, start_all, units='mins')),
                                as.numeric(difftime(end_fit, start_fit, units='mins')),
                                as.numeric(difftime(end_est, start_est, units='mins'))
                    )  )
  
  res <- est_df %>% 
    full_join(diags_df, by=c("estimand","value")) %>%
    full_join(timings, by=c("estimand","value") )
  
  print(tail(res)) 
  
  if (!is.null(save_file)) saveRDS(res, file=save_file)
  
  
}





comb_res_slurm <- function(i){
  
  est <- NULL
  
  for(j in 1:nsim){
    
    if(file.exists(paste0("../scen", i, "/est", j,".rds"))){
      temp <- readRDS(paste0("../scen", i, "/est", j,".rds"))
      temp <- temp %>%
        mutate(isim = j, mod_type = "sim") 
      est <- rbind(est, temp)
    }
    
    
  }
  
  saveRDS(est, paste0("../scen", i, "_est.rds"))
  
  true <- mod_est(mod = true_mod, estimand = NULL, t_vec, rmst_vec, nonprop = FALSE,
                  return_model = TRUE, save_file = NULL)
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("../scen", i, "_res.rds"))
  
}



## Evaluate estimands for true/fitted models

#' @param mod survextrap/flexsurv object containing fitted model
#' @param path_to_store store data in "data" directory
#' @param estimands list of estimands to evaluate.
#' @param t, time vector for time-varying estimands such as e.g. hazard function
#' @param rmst_t, tvec for rmst calculation.
#

mod_est <- function(mod, estimands = NULL, t, rmst_t, nonprop = FALSE,
                    return_model = FALSE, save_file = NULL){
  
  model_class <- class(mod)

  if(model_class == "survextrap"){
   
  summ_fns <- list(median=median, ~quantile(.x, probs=c(0.025, 0.975)), se = sd)
    
  rmst_survextrap <- rmst(mod, t=rmst_t,  summ_fns=summ_fns) 
  haz_survextrap <- hazard(mod, t=t, summ_fns=summ_fns)
  surv_survextrap <- survival(mod, t=t, summ_fns=summ_fns)

  res <- 
    tibble(
      estimand = "rmst",
      t = rmst_t,
      value = rmst_survextrap[["median"]],
      value_ci_low = rmst_survextrap[["2.5%"]],
      value_ci_high = rmst_survextrap[["97.5%"]],
      value_se = rmst_survextrap[["se"]]
    )      %>%
    bind_rows(
      tibble(
        estimand = "hazard",
        t = t,
        value = haz_survextrap[["median"]],
        value_ci_low = haz_survextrap[["2.5%"]],
        value_ci_high = haz_survextrap[["97.5%"]],
        value_se = haz_survextrap[["se"]]
      )) %>%
    bind_rows(
      tibble(
        estimand = "survival",
        t = t,
        value = surv_survextrap[["median"]],
        value_ci_low = surv_survextrap[["2.5%"]],
        value_ci_high = surv_survextrap[["97.5%"]],
        value_se = surv_survextrap[["se"]] 
      ))          

  }
  
  if(model_class == "flexsurvreg"){
    
  
  rmst_flexsurvreg <- standsurv(mod, t=rmst_t, type = "rmst", ci = TRUE, se = TRUE, boot = T, B = 1000)
  haz_flexsurvreg <- standsurv(mod, t=t, type = "hazard", ci = TRUE, se = TRUE, boot = T, B = 1000)
  surv_flexsurvreg <- standsurv(mod, t=t, type = "survival", ci = TRUE, se = TRUE, boot = T, B = 1000)
  
  res <- 
    tibble(
      estimand = "rmst",
      t = rmst_t,
      value = rmst_flexsurvreg[["at1"]],
      value_ci_low = rmst_flexsurvreg[["at1_lci"]],
      value_ci_high = rmst_flexsurvreg[["at1_uci"]],
      value_se = rmst_flexsurvreg[["at1_se"]]
    )      %>%
    bind_rows(
      tibble(
        estimand = "hazard",
        t = t,
        value = haz_flexsurvreg[["at1"]],
        value_ci_low = haz_flexsurvreg[["at1_lci"]],
        value_ci_high = haz_flexsurvreg[["at1_uci"]],
        value_se = haz_flexsurvreg[["at1_se"]]
      )) %>%
    bind_rows(
      tibble(
        estimand = "survival",
        t = t,
        value = surv_flexsurvreg[["at1"]],
        value_ci_low = surv_flexsurvreg[["at1_lci"]],
        value_ci_high = surv_flexsurvreg[["at1_uci"]],
        value_se = surv_flexsurvreg[["at1_se"]]
      ))          
  
  }
  
  if(!is.null(save_file)) saveRDS(res, file = save_file)
  
  if(return_model == TRUE)  res
  
}

##### treatment version

mod_est_trt <- function(mod, estimands = NULL, t, rmst_t, nonprop = FALSE,
                    return_model = FALSE, save_file = NULL){
  
  model_class <- class(mod)
  
  if(model_class == "survextrap"){
    
    summ_fns <- list(median=median, ~quantile(.x, probs=c(0.025, 0.975), na.rm = TRUE), se = sd)
    
    newdata_trt <- tibble(trt = as.factor(c(0,1)))
    
    rmst_survextrap <- rmst(mod, t=rmst_t, newdata = newdata_trt, summ_fns=summ_fns)  %>%
      mutate(estimand = "rmst") %>%
      select(estimand, trt, t, median, `2.5%`,  `97.5%`, se ) %>%
      rename(value = median, 
             value_ci_low = `2.5%`,
             value_ci_high =`97.5%`,
             value_se = se)
    
    irmst_survextrap <- irmst(mod, t=rmst_t, newdata = newdata_trt, summ_fns=summ_fns)  %>%
      mutate(estimand = "irmst", trt = NA) %>%
      select(estimand, trt, t, median, `2.5%`,  `97.5%`, se ) %>%
      rename(value = median, 
             value_ci_low = `2.5%`,
             value_ci_high =`97.5%`,
             value_se = se)

    hr_survextrap <- hazard_ratio(mod, t=t, newdata = newdata_trt, summ_fns=summ_fns) %>%
      mutate(estimand = "hr", trt = NA) %>%
      select(estimand, trt, t, median, `2.5%`,  `97.5%`, se ) %>%
      rename(value = median, 
             value_ci_low = `2.5%`,
             value_ci_high =`97.5%`,
             value_se = se)
    
    haz_survextrap <- hazard(mod, t=t, newdata = newdata_trt, summ_fns=summ_fns) %>%
      mutate(estimand = "hazard") %>%
      select(estimand, trt, t, median, `2.5%`,  `97.5%` , se) %>%
      rename(value = median, 
             value_ci_low = `2.5%`,
             value_ci_high =`97.5%`,
             value_se = se)

    surv_survextrap <- survival(mod, t=t, newdata = newdata_trt, summ_fns=summ_fns) %>%
    mutate(estimand = "survival") %>%
      select(estimand, trt, t, median, `2.5%`,  `97.5%` , se) %>%
      rename(value = median, 
             value_ci_low = `2.5%`,
             value_ci_high =`97.5%`,
             value_se = se)

    res <- rmst_survextrap   %>%
      bind_rows(irmst_survextrap) %>%
      bind_rows(hr_survextrap) %>%
      bind_rows(haz_survextrap) %>%
      bind_rows(surv_survextrap)

  }
  
  if(model_class == "flexsurvreg"){
    
    
   print("trt needs implementing...")
    
  }
  
  
  
  if(!is.null(save_file)) saveRDS(res, file = save_file)
  
  if(return_model == TRUE)  res
  
}

# Function to combine together results across simulations. 

comb_res_slurm <- function(i){
  
  est <- NULL
  
  for(j in 1:nsim){
    
    if(file.exists(paste0("../scen", i, "/est", j,".rds"))){
      temp <- readRDS(paste0("../scen", i, "/est", j,".rds"))
      temp <- temp %>%
        mutate(isim = j, mod_type = "sim") 
      est <- rbind(est, temp)
    }
    
    
  }
  
  saveRDS(est, paste0("../scen", i, "_est.rds"))
  
  true <- mod_est(mod = true_mod, estimand = NULL, t_vec, rmst_vec, nonprop = FALSE,
                  return_model = TRUE, save_file = NULL)
  
  true <- true %>%
    mutate(isim = 0, mod_type = "true")
  
  res <- create_res(true, est)
  
  saveRDS(res, paste0("../scen", i, "_res.rds"))
  
}


# function to create simulation results tibble 
# to feed into performance.R
# also add IDs to the estimands.
# 

create_res <- function(true, est){
  
  res <- rbind(true, est) %>%
      group_by(estimand, isim) %>%
      mutate(id = row_number(),
             estimand_id = paste0(estimand, id)) %>%
      relocate(estimand_id) %>%
      select(-id) %>%
      ungroup()
  
  res
  
}

# Create ....

diags <- function(mod){
  if (!(mod$fit_method=="mcmc")){
    dg <- list(prop_divergent=NA, warning_any=NA, warning_divergent=NA,
               warning_treedepth=NA, warning_BayesianFrac=NA,
               warning_SamplingProb=NA, warning_rhatHigh=NA,
               warning_BulkEssLow=NA, warning_TailEssLow=NA)
  } else {
    dg <- throw_sampler_warnings(mod$stanfit) 
    dg$warning_any <- 1*(sum(dg != 0) > 0)
    
    chain_warnings <- rstan::get_sampler_params(mod$stanfit, inc_warmup=FALSE)
    div <- NULL
    for(i in 1:length(chain_warnings)){
      div_chain <- chain_warnings[[i]][,'divergent__']
      div <- c(div, div_chain)
    }
    dg$prop_divergent <- sum(div)/length(div)
  }
  dg <- tibble::enframe(dg) |>
    tidyr::unnest(cols="value")
  dg
}


# test model divergence:
catchToList <- function(expr) {
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
  list(value = val, warnings = myWarnings, error=myError)
} 


# stanfit object

throw_sampler_warnings <- function(object) {
  
  divergent <- 0
  treedepth <- 0
  BayesianFrac <- 0
  SamplingProb <- 0
  rhatHigh <- 0
  BulkEssLow <- 0
  TailEssLow <- 0
  
  if (!is(object, "stanfit"))
    stop("'object' must be of class 'stanfit'")
  
  if (isTRUE(object@stan_args[[1]][["algorithm"]] == "Fixed_param")) {
    return(invisible(NULL))
  }
  
  sp <- get_sampler_params(object, inc_warmup = FALSE)
  n_d <- sum(sapply(sp, FUN = function(x) {
    if ("divergent__" %in% colnames(x)) return(sum(x[,"divergent__"]))
    else return(0)
  }))
  if (n_d > 0) {
    ad <- object@stan_args[[1]]$control
    if (is.null(ad)) ad <- 0.8
    else {
      ad <- ad$adapt_delta
      if (is.null(ad)) ad <- 0.8
    }
    divergent <- 1
  }
  max_td <- object@stan_args[[1]]$control
  if (is.null(max_td)) max_td <- 10
  else {
    max_td <- max_td$max_treedepth
    if (is.null(max_td)) max_td <- 10
  }
  n_m <- sum(sapply(sp, FUN = function(x) {
    if ("treedepth__" %in% colnames(x)) return(sum(x[,"treedepth__"] >= max_td))
    else return(0)
  }))
  if (n_m > 0)
    treedepth <- 1
  n_e <- 0L
  if (all(sapply(sp, function(x) "energy__" %in% colnames(x)))) {
    E <- as.matrix(sapply(sp, FUN = function(x) x[,"energy__"]))
    threshold <- 0.2
    if (nrow(E) > 1) {
      EBFMI <- get_num_upars(object) / apply(E, 2, var)
      n_e <- sum(EBFMI < threshold, na.rm = TRUE)
    }
    else n_e <- 0L
    if (n_e > 0)
      BayesianFrac <- 1
  }
  if (n_d > 0 || n_m > 0 || n_e > 0) 
    SamplingProb <- 1
  
  sims <- as.array(object)
  rhat <- apply(sims, MARGIN = 3, FUN = Rhat)
  if (any(rhat > 1.05, na.rm = TRUE))
    rhatHigh <- 1
  bulk_ess <- apply(sims, MARGIN = 3, FUN = ess_bulk)
  if (any(bulk_ess < 100 * ncol(sims), na.rm = TRUE))
    BulkEssLow <- 1
  tail_ess <- apply(sims, MARGIN = 3, FUN = ess_tail)
  if (any(tail_ess < 100 * ncol(sims), na.rm = TRUE))
    TailEssLow <- 1
  
  return(list(divergent=divergent,
              treedepth=treedepth,
              BayesianFrac=BayesianFrac,
              SamplingProb=SamplingProb,
              rhatHigh=rhatHigh,
              BulkEssLow=BulkEssLow,
              TailEssLow=TailEssLow))
}

