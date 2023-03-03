#Subject: survival hazard function in R
#https://devinincerti.com/code/survival-distributions.html

#Use “hgamma” function in package “flexsurv”

install.packages("flexsurv")
library(flexsurv)



#define gamma hazard fn: shape parameter (gamma), scale(rate) parameter (lambda), rate=1 for standard gamma distribution

gamma_haz <- function(t, x, betas, shape, rate = 1){
  
  exp(betas * x) * flexsurv::hgamma(t, shape = shape, rate = rate)
  
}

#generate mixture data

covs <- data.frame(id = 1:N,
                   trt = stats::rbinom(N, 1, 0.5))

dat <- simsurv(hazard = gamma_haz,
               shape = 2,
               rate = 1,
               betas = c(trt = 2),
               x = covs,
               maxt = 5)


library(tydiverse)
# modify the gamma hazard function to force it returns values of opposite sign at the endpoints to ensure that the uniroot function can find a root
gamma_haz2 <- function(t, x, betas, shape, rate = 1) {
  h <- exp(betas * x) * flexsurv::hgamma(t, shape = shape, rate = rate)
  if (length(h) > 1) {
    h[1] <- -h[1]
    h[length(h)] <- h[length(h)]
  } else {
    h <- ifelse(t > 0, -h, h)
  }
  return(h)
}
# gamma hazard fn return to a single value for each time point
gamma_haz3 <- function(t, x, betas, shape, rate = 1) {
  h <- exp(betas * x) * flexsurv::hgamma(t, shape = shape, rate = rate)
  ifelse(t == 0, -h, h)
}

#write a fn to simulate gamma data
sim_gamma <- function(k, N){
  #generate exponential data
  covs <- data.frame(id = 1:N,
                     trt = stats::rbinom(N, 1, 0.5))
  dat <- cbind(simsurv(hazard = gamma_haz2,
                       shape = 2,
                       rate = 1,
                       betas = c(trt = 2),
                       x = covs,
                       maxt = 5), covs)
  
  #fit models
  fit.exponential <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "exponential")
  fit.weibull <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "weibull")
  fit.cox <- coxph(Surv(eventtime, status) ~ trt, data = dat)
  
  #extract beta
  result <- tibble(exp_beta = c(-fit.exponential$coefficients[-1]), 
                   weibull_beta = c(-fit.weibull$coefficients[-1])/fit.weibull$scale,
                   cox_beta = c(fit.cox$coefficients), 
                   dist = "gamma",
                   beta = 2, 
                   gamma = 2,
                   N = N)
  return(result)
}

# Set seed for reproducibility
set.seed(2023)

# Create empty dataframe to store results
sim_gamma_result  <- data.frame()

# Simulate 1000 times
# for (n in c(100, 150, 200, 250, 300, 350, 400)) {
# for (i in 1:1000)

for (n in c(100,150)) {
  for (i in 1:10) {
    sim_res <- sim_gamma(N = n)
    sim_gamma_result <- rbind(sim_gamma_result, sim_res)
  }
}

gamma_table <- sim_gamma_result %>% 
  group_by(N) %>%
  summarize(mse_exp = mean((exp_beta-2)^2),
            mse_weibull = mean((weibull_beta-2)^2),
            mse_cox = mean((cox_beta-2)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta-2),
            bias_weibull = mean(weibull_beta-2),
            bias_cox = mean(cox_beta-2)
            
  )

gamma_table

gamma_table <- gamma_table %>% 
  pivot_longer(mse_exp:bias_cox, values_to = "performance", names_to = "method") %>% separate_wider_delim(method, "_", names = c("method", "model")) %>%
  pivot_wider(names_from = method, values_from = performance)
gamma_table
write.csv(gamma_table,"table/gamma_result_table.csv")


