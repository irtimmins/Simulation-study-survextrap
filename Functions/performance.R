## evaluate estimands survextrap model

#' @param res parameter estimates/output for true and simulated data
#' @param estimands vector of estimands to be assessed

## for review

sim_perform <- function(res, estimands = c("rmst")){
  
  # extract results for estimands of interest
  res_par <- res %>%
    filter(estimand %in% estimands) 
  
  # extract simulation replicates
  sim_par <- res_par %>%
    filter(mod_type == "sim") 

  # extract true parameter values for reference
  true_par <- res_par %>%
    filter(mod_type == "true") 
  
  # empty tibble to store results
  per_summ_all <- NULL
  
  for(i in 1:nrow(true_par)){
        
        # extract estimand id (labelled in create_res function)
        est_id <- true_par$estimand_id[i]
        
        # true parameter value for assessing bias
        true_val <- true_par$value[i]

        # extract simulation replicates for input into simsum
        simsum_in <- sim_par %>%
          filter(estimand_id == est_id)
        
        per_out <- simsum(data = simsum_in, estvarname = "value",
                          true = true_val, 
                          se = "value_se", ci.limits = c("value_ci_low", "value_ci_high"))
        per_summ <- per_out$summ %>% mutate(estimand_id = est_id)
        
        if (i == 1) per_summ_all <-  per_summ
        else per_summ_all <-  bind_rows(per_summ_all, per_summ)

  }

  as_tibble(per_summ_all)
}



