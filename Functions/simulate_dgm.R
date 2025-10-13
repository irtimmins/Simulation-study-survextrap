## Simulate from the "true" survival model

#' @param model_class class of true model e.g. flexsurvreg
#' @param true_mod true model 
#' @param N Number of patients in trial
#' @param nsim number of simulated datasets
#' @param maxT maximum event time
#' @param seed Random number seed

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


# Simulate data for two arm DGMs.

#' @param model_class class of true model e.g. flexsurvreg
#' @param true_mod true model 
#' @param N Number of patients in trial
#' @param nsim number of simulated datasets
#' @param maxT maximum event time
#' @param seed Random number seed

sim_dgm_trt <- function(gamma = NULL, 
                        knots = NULL,
                        hr_scenario = NULL,
                        beta = NULL,
                        N,
                        nsim, 
                        maxT,
                        seed,
                        lower = 1e-08, 
                        upper = 10000, 
                        nodes=100){
  
  if(hr_scenario == 1) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      t <- as.matrix(t)
      trow <- nrow(t)
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+beta[1]*trt
      loghaz <- matrix(loghaz, nrow = trow)
      return(exp(loghaz))
    }
    
  } else if(hr_scenario == 2) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      t <- as.matrix(t)
      trow <- nrow(t)
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+
        beta[1]*(1-tanh(x = (beta[2]*t+beta[3])))*trt
      loghaz <- matrix(loghaz, nrow = trow)
      return(exp(loghaz))
    } 
    
  } else if(hr_scenario == 3) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      t <- as.matrix(t)
      trow <- nrow(t)
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+
        (beta[1]+beta[2]*demg(x = t, mu = beta[3], sigma = beta[4], lambda = beta[5]))*trt
      loghaz <- matrix(loghaz, nrow = trow)
      return(exp(loghaz))
    }  
    
  } else if(hr_scenario == 4) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      
      s <- function(t) 0.5*(1 + tanh(beta[1]*(t - 1)))
      d <- function(t) 1 - 1/(1 + exp(-beta[2]*(t - 2.5)))
      
      t <- as.matrix(t)
      trow <- nrow(t)
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+ 
        (beta[3] + (beta[4] - beta[3]) * s(t)) * d(t)*trt
      
      loghaz <- matrix(loghaz, nrow = trow)
      return(exp(loghaz))
    }  
  
  
  ## Cumulative hazard function using Guassian quadrature.
  
  f1 <- function(t,lnU,nodes,gamma, knots, beta, trt) {
    ch <- vecquad_gl(haz,0,t,nodes, gamma, knots, beta, trt)
    return(-ch-lnU)
  }
  
  vecquad_gl <- function(fn,a,b,nodes,...) {
    nw  <- gaussLegendre(n = nodes,-1,1)
    z1 <- fn(t = 0.5*as.matrix(b-a)%*%t(nw$x)+as.vector(0.5*(a+b)),...)
    
    return((0.5*(b-a))*rowSums(sweep(z1,2, nw$w,"*")))
  }  
  
  
  set.seed(seed)
  
  sim_data <- expand_grid(id=1:floor(N/2),
                          censtime = maxT,
                          i=1:nsim,
                          trt = c(0,1)) %>%
    mutate(time = 0,
           event = 0)
  
  root_interval <- cbind(rep(lower,nrow(sim_data)),rep(upper,nrow(sim_data)))
  lnU  <- log(runif(nrow(sim_data)))
  
  t <- vuniroot(f1,interval=root_interval, lnU=lnU, nodes=nodes, 
                gamma = gamma, knots = knots,  beta = beta, trt = sim_data[["trt"]])$root
  
  
  if(is.null(maxT)) maxT <- Inf
  
  sim_data <- sim_data %>%
    mutate(time = pmin(t, censtime),
           event = as.numeric(t<censtime))
  
  return(sim_data)
  
}


