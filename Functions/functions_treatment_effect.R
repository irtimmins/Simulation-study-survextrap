
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
    
  }
  
}
