## evaluate estimands for true/fitted models

#' @param mod survextrap/flexsurv object containing fitted model
#' @param path_to_store store data in "data" directory
#' @param estimands list of estimands to evaluate.
#' @param t, time vector for time-varying estimands such as e.g. hazard function
#' @param rmst_t, tvec for rmst calculation.
#

## for review

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

  ## add HRs fo treatment effects...
  
  #if(nonprop == FALSE){
  
  #loghr_est =  # extract from survextrap summary   
  
  #}else{
  
  # hr_nonprop_est = # Hazard Ratios for each time in tvec
  
  #}
  
  ## add stan algorithm performance characteristics...
  
  # stan_extract <- function(fit_mod){...} 
  # stan_est <- stan_extract(fit_mod)
  
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
      
    #print(head(res))
    
    ## add HRs fo treatment effects...
    
    #if(nonprop == FALSE){
    
    #loghr_est =  # extract from survextrap summary   
    
    #}else{
    
    # hr_nonprop_est = # Hazard Ratios for each time in tvec
    
    #}
    
    ## add stan algorithm performance characteristics...
    
    # stan_extract <- function(fit_mod){...} 
    # stan_est <- stan_extract(fit_mod)
    
  }
  
  if(model_class == "flexsurvreg"){
    
    
   print("trt needs implementing...")
    
  }
  
  
  
  if(!is.null(save_file)) saveRDS(res, file = save_file)
  
  if(return_model == TRUE)  res
  
}






# for review - although combine_est may be replaced more by slurm parallelisation
# combine estimates across simulations.

combine_est <- function(path_to_mod, nsim = nsim,
                        estimands = NULL, t, rmst_t, nonprop = FALSE,
                        return_model = FALSE, save_file = NULL){

comb <- NULL  
  
  for(i in 1:nsim){

    file_mod <- paste0(path_to_mod, i,".rds")
    file_store <- paste0(save_file, i,".rds")
    
    mod <- readRDS(file_mod)

    temp <- mod_est(mod, estimands=NULL, t, rmst_t, nonprop=NULL,
                    return_model=TRUE, save_file = file_store)
    
    temp <- temp %>%
      mutate(isim = i, mod_type = "sim")
    
    comb <- rbind(comb, temp)
    
    
  }

  comb

}

# for review

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








#' @param t_chr, time vector for time-varying estimands such as e.g. hazard function
#' @param rmst_t_chr, tvec for rmst calculation.
#

sim_est_slurm <- function(fit_mod_path, estimand = NULL, t_chr, rmst_t_chr, nonprop = FALSE,
                          return_model = FALSE, save_file = NULL ){
  
  
  rmst_t <- get(rmst_t_chr)
  t <- get(t_chr)
  mod <- readRDS(fit_mod_path)
  #print(summary(mod))
  mod_est(mod, estimand = NULL, t, rmst_t, nonprop = FALSE,
          return_model, save_file )
  
}

