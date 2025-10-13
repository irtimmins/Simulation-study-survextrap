
#######################################################
# Hazard-ratio and Hazard functions for nonPH scenarios.
#######################################################


scenario_hr <- function(hr_scenario){
  
  if(hr_scenario == 1) {
    
    hr <- function(t, gamma, knots, beta) {
      loghr <- rep(beta[1], length(t))
      return(exp(loghr))
    }
    
  } else if(hr_scenario == 2) {
    
    hr <- function(t, gamma, knots, beta) {
      loghr <- beta[1]*(1-tanh(x = (beta[2]*t+beta[3])))
      return(exp(loghr))
    } 
    
  } else if(hr_scenario == 3) {
    
    hr <- function(t, gamma, knots, beta) {
      loghr <- (beta[1]+beta[2]*demg(x = t, mu = beta[3], sigma = beta[4], lambda = beta[5]))
      return(exp(loghr))
    }  
    
  } else if(hr_scenario == 4) {
    
    hr <- function(t, gamma, knots, beta) {
      
      s <- function(t) 0.5*(1 + tanh(beta[1]*(t - 1)))
      d <- function(t) 1 - 1/(1 + exp(-beta[2]*(t - 2.5)))
      
      loghr <- (beta[3] + (beta[4] - beta[3]) * s(t)) * d(t)
      
      return(exp(loghr))
    }  
    
  }
}
  
  
scenario_haz <- function(hr_scenario){
  
  if(hr_scenario == 1) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+
        beta[1]*trt
      return(exp(loghaz))
    }
    
  } else if(hr_scenario == 2) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+ 
        beta[1]*(1-tanh(x = (beta[2]*t+beta[3])))*trt
      return(exp(loghaz))
    } 
    
  } else if(hr_scenario == 3) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+
        (beta[1]+beta[2]*demg(x = t, mu = beta[3], sigma = beta[4], lambda = beta[5]))*trt
      return(exp(loghaz))
    }  
    
  }else if(hr_scenario == 4) {
    
    haz <- function(t, gamma, knots, beta, trt) {
      
      s <- function(t) 0.5*(1 + tanh(beta[1]*(t - 1)))
      d <- function(t) 1 - 1/(1 + exp(-beta[2]*(t - 2.5)))
      
      loghaz <- log(hsurvspline(x = t, gamma = gamma, knots = knots))+ 
        (beta[3] + (beta[4] - beta[3]) * s(t)) * d(t)*trt

      return(exp(loghaz))
    }  
    
  }
  
}
