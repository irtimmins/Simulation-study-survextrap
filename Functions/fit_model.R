# Fit survextrap model to each simulation replicate
#'
#' @param sim_data character name of simulated individual level data, from sim_dgm()
#' @param prior character name of prior object
#' @param mspline character name of mspline object
#' @param stan_fit character name of stan parameters object
#' @param maxT maximum time for hazard rate calculation
#' @param return_model  if \code{return_model=TRUE} then return object in local environment/memory 
#' @param save_file then write external copy of fitted model to specified file location
#' 

# for review - though fit_model_all may be replaced by slurm parallelisation

fit_model_all <- function(sim_data, prior, mspline, stan_fit, maxT,
                          return_model = FALSE, save_file = NULL,
                          seed = NULL){
  
  nsim <- max(sim_data$i)
  
  for(j in 1:nsim){

    save_file_j <- NULL
    
    if(!is.null(save_file))  save_file_j <- paste0(save_file, j,".rds")
    
    temp <- extract_data(sim_data = sim_data, 
                         isim = j)
    
    fit_temp <- fit_model(temp,
                          prior,
                          mspline, 
                          stan_fit, 
                          maxT,
                          return_model,
                          save_file = save_file_j,
                          seed)
    
  }
  
}


#' @param data simulated individual level data, for single trial data.
#' @param prior list of survextrap prior parameters
#' @param mspline list of mspline parameters
#' @param stan_fit  list of stan algorithm parameters 

# for review

fit_model <- function(data, prior, mspline, stan_fit, maxT, 
                      return_model = FALSE, save_file = NULL ,
                      seed = NULL){
  
  set.seed(seed)
  
  res <- survextrap(Surv(time, event) ~ 1, 
                    mspline = mspline, 
                    prior_hscale = prior$prior_hscale,
                    prior_hsd = prior$prior_hsd,
                    data = data,
                    fit_method=stan_fit$fit_method)           

  if(!is.null(save_file)) saveRDS(res, file = save_file)
  if(return_model == TRUE)  res
  
}


extract_data <-  function(sim_data, isim){
  data <-  sim_data %>%
    filter(i == isim) %>%
    select(-i)
  data
}


#' @param sim_data_chr character name of simulated individual level data, from sim_dgm()
#' @param prior_chr character name of prior object
#' @param mspline_chr character name of mspline object
#' @param stan_fit_chr character name of stan parameters object

## (not for review) fit_model_slurm function is still under development

fit_model_slurm <- function(sim_data_chr, isim, prior_chr, mspline_chr, stan_fit_chr, maxT,
                            return_model = FALSE, save_file = NULL, seed = NULL){
  
  sim_data <- get(sim_data_chr)
  prior <- get(prior_chr)
  mspline <- get(mspline_chr)
  stan_fit <- get(stan_fit_chr)

  data <- extract_data(sim_data, isim)

  
  fit_model(data, prior, mspline, stan_fit, maxT,  
            return_model, save_file, seed)
  
  
}


hazard_area_diff <- function(x, maxt=NULL, mint=0, haz_ref, diff_type=2, ...){
  summ_mod <- summary(x)
  
  mod_knots <- x$mspline$knots
  
  mod_alpha <- summ_mod %>% 
    filter(variable == "alpha") %>%
    select(median) %>% pull()
  
  mod_coefs <- summ_mod %>% 
    filter(!(variable %in% c("alpha", "hsd"))) %>%
    select(median) %>% pull()
  
  fn <- function(t){
    haz <- hsurvmspline(x = t, 
                        alpha = mod_alpha,
                        coefs = mod_coefs,
                        knots = mod_knots)
    
    if (diff_type == 1)
      ret <- abs(haz - haz_ref(t))
    else if (diff_type == 2)
      ret <- (haz - haz_ref(t))^2
    else stop("diff_type should be 1 or 2")
  }
  int <- try(integrate(fn, lower=mint, upper = maxt, subdivisions = 5000))
  
  int
}


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












