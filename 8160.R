library(survival)
library(tidyverse)
library(fastmap)
set.seed(2023)
simdata = function(n, distribution,b1,lambda=0.5, gamma=2,alpha=0.5) {
  exponential = function(u, x, lambda, b1) {
    time = -log(u) / (exp(x * b1) * lambda) }  
  weibull = function(u, x, lambda, gamma, b1) {
    time = ( -log(u) / (exp(x * b1) * lambda) ) ^ (1 / gamma) }
  gompertz = function(u, x, lambda, alpha, b1) {
    time = (1/alpha) * log( 1 - (alpha * log(u) / (lambda * exp(b1 * x))) ) }
  x = rbinom(n, 1, 0.5)
  u = runif(n)

  if (distribution == "exponential") {
    time = exponential(u, x, lambda, b1) } 
  else if (distribution == "weibull") {
    time = weibull(u, x, lambda, gamma,b1) }
  else { 
    time = gompertz(u, x, lambda, alpha, b1) } 
  
  time[time < 1/365] = 1/365
  time[time == 1 / 365] = time[time == 1 / 365] + rnorm(length(time[time == 1 / 365]), 0, 1e-4)
  time = abs(time)
  e=as.numeric(time<5)
  time=pmin(time,5)
  
  # data set
  survival_data = data.frame(id=1:n,
                             time = time,
                             x = x,
                             event=e)
}


set.seed(2023)
simulate = function(sim, n, b1, dist = "weibull") {
  # Set up coefficient vectors
  exp_beta = rep(NA, sim)
  weibull_beta = rep(NA, sim)
  gompertz_beta = rep(NA, sim)
  
 
    
    for (i in 1:sim) {
      
      data = simdata(n, distribution = dist,b1, 
                           lambda = 0.5, gamma = 2)
      
      fit.exponential = survreg(Surv(data$t, data$event) ~ data$x, 
                                dist = "exponential") 
      fit.weibull = survreg(Surv(data$t, data$event) ~ data$x, 
                            dist = "weibull")
      fit.gompertz = coxph(Surv(data$t, data$event) ~ data$x)
      
       
      exp_beta[i] = -fit.exponential$coefficients[-1]
      weibull_beta[i] = -fit.weibull$coefficients[-1] / fit.weibull$scale
      gompertz_beta[i] = fit.gompertz$coefficients[1]
      
    }
    
   
  
  coef = tibble(exp = exp_beta,
                weibull = weibull_beta,
                gompertz=gompertz_beta)
  exp_var = var(exp_beta)
  weibull_var = var(weibull_beta)
  gompertz_var=var(gompertz_beta)
  
  char = tibble(model = c("exp", "weibull","gompertz"),
                var=c(exp_var,weibull_var,gompertz_var))
  
  results = list(coef, char)
  names(results) = c("coefficients", "performance")
  #return(results)
  return(coef)
  
  }
  

set.seed(2023)
sim_results = 
  tibble(sample_size = c(100,150,200,250,300,350,400)) %>% 
  mutate(
    output= map(.x = sample_size, ~simulate(sim = 1000, n = .x, b1 = -0.5, 
                                                        dist = "weibull")$performance),
                output = map(.x = output, ~mutate(.x)))

fit_exp = function(df) {
  fit.exponential = survreg(Surv(time, event) ~ x, dist = "exponential", data = df)
  return(as.numeric(-fit.exponential$coefficients[-1]))
}

fit_weibull = function(df) {
  fit.weibull <- survreg(Surv(time, event) ~ x, dist = "weibull", data = df)
  return(as.numeric(-fit.weibull$coefficients[-1] / fit.weibull$scale))
}

fit_cox = function(df) {
  fit.cox <- coxph(Surv(time, event) ~ x, data = df)
  return(as.numeric(fit.cox$coefficients))
}

set.seed(2023)
sim_beta=expand.grid(sample_size=c(100,150,200,250,300,350,400),rep=1:1000) %>%
  mutate(data=map(.x=sample_size,~simdata(n = .x,distribution="weibull", b1=-0.5))) %>%
  mutate(
    exp_beta = map_dbl(data, fit_exp),
    weibull_beta = map_dbl(data, fit_weibull),
    cox_beta = map_dbl(data, fit_cox)
  ) 

set.seed(2023)
sim_beta %>%ggplot()+geom_density(aes(x = exp_beta,col="red"), size = 0.7) +
  geom_density(aes(x = weibull_beta,col="blue"), size = 0.7) +
  geom_density(aes(x = cox_beta,col="yellow"),size = 0.7)+
  facet_grid(sample_size~.)+scale_color_manual(name="models",values=c("red"="red","yellow"="yellow","blue"="blue"),
                                               labels=c("gompertz","exp","weibull"))+
  scale_x_continuous(breaks = c(-1, -0.5, 0), limits = c(-1.5, 0.5))


temp1 = sim_results %>% unnest(output) 

ggplot(temp1, aes(x = sample_size, y = var, color = model)) +
  geom_point(size = 1) +
  geom_line(size = 1) +
  xlab("sample size") +
  labs(caption = "Figure 1",
       title="Data from Weibull distribution fits in 3 models") + 
  theme(plot.caption = element_text(hjust = 0.5, size = rel(1.2)))


#confidence interval
set.seed(2023)
g1=simulate(sim=1000,n=400,b1=-0.5)
g1
exp1=t.test(g1$exp,conf.level = 0.95)
weibull1=t.test(g1$weibull,conf.level = 0.95)
weibull1
cox1=t.test(g1$gompertz,conf.level = 0.95)
cox1
