library(survival)
library(tidyverse)

set.seed(2023)
simdata = function(n, distribution,b1,lambda=0.5, gamma=2) {
  exponential = function(u, x, lambda, b1) {
    time = -log(u) / (exp(x * b1) * lambda) }  
  weibull = function(u, x, lambda, gamma, b1) {
    time = ( -log(u) / (exp(x * b1) * lambda) ) ^ (1 / gamma) }
  x = rbinom(n, 1, 0.5)
  u = runif(n)

  if (distribution == "exponential") {
    time = exponential(u, x, lambda, b1) } 
  else{
   time = weibull(u, x, lambda, gamma, b1)} 
  
  e = runif(n)
  t= pmin(time, e)
  event=ifelse(time<=e,1,0)
  
  # data set
  survival_data = data.frame(id=1:n,
                             time = time,
                             t = t,
                             x = x,
                             event=event)
}


fit_exp = function(df) {
  fit.exponential = survreg(Surv(t, event) ~ x, dist = "exponential", data = df)
  return(as.numeric(-fit.exponential$coefficients[-1]))
}



fit_weibull = function(df) {
  fit.weibull <- survreg(Surv(t, event) ~ x, dist = "weibull", data = df)
  return(as.numeric(-fit.weibull$coefficients[-1] / fit.weibull$scale))
}


simulate = function(sim, n, b1, dist = "exponential") {
  # Set up coefficient vectors
  exp_beta = rep(NA, sim)
  weibull_beta = rep(NA, sim)
  cox_beta = rep(NA, sim)
  
 
    
    for (i in 1:sim) {
      # Generate data
      data = simdata(n, distribution = dist,b1, 
                           lambda = 0.5, gamma = 2)
      # Fit three survival distributions
      fit.exponential = survreg(Surv(data$t, data$event) ~ data$x, 
                                dist = "exponential") 
      fit.weibull = survreg(Surv(data$t, data$event) ~ data$x, 
                            dist = "weibull")
      
      # Save beta coefficients 
      exp_beta[i] = -fit.exponential$coefficients[-1]
      weibull_beta[i] = -fit.weibull$coefficients[-1] / fit.weibull$scale
      
    }
    
   
  # Store beta coefficients
  coef = tibble(exp = exp_beta,
                weibull = weibull_beta,
                ) }
  
  si1=simulate(500,500,0.5)


