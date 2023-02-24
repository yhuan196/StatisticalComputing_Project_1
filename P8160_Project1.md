P8160_Project1_Survival
================
Yi Huang, yh3554
2023-02-19

# Project 1: Design a simulation study to compare three survival models

## Examples : Simsurv Package Data Generarion

``` r
set.seed(2023)
# tot number of patients
N <- 1000

# define covariates
covs <- data.frame(id = 1:N, trt = stats::rbinom(N, 1, 0.5))

## Exponential
exp_dist <- simsurv(dist = "exponential", lambdas = 0.5, 
                    x = covs, betas = c(trt = -0.5), maxt = 5)


## Weibull
weibull_dist <- simsurv(dist <- "weibull", lambdas = 0.5, gammas = 0.05, 
                        x = covs, betas = c(trt = -0.5), maxt = 5)


## Gompertz
gompertz_dist <- simsurv(dist = "gompertz", lambdas = 0.1, gammas = 0.05, 
                         x = covs, betas = c(trt = -0.5), maxt = 5)
```

## Generate Gompertz distribution use inverse transformation method

``` r
gen_gompertz <- function(alpha = 0.5, lambda = 0.5, b1 = -0.5, n, seed) {
  # Generate a random sample from the uniform distribution
  set.seed(seed)
  u <- runif(n)
  x <- rep(0, n)
  x <- (1/alpha)*log(1-(alpha*log(u)/(lambda*exp(x*b1)))) 

  return(x)
  
}

sample_gompertz <- gen_gompertz(n=1000, seed = 123123)
df <- data.frame(sample_gompertz)

# Visualization to validate the algorithm
 df %>% ggplot(aes(x = sample_gompertz)) +
  geom_histogram(aes(y = after_stat(density)),
                 colour = 1, fill = "red", bins = 30) +
  geom_density() + xlab("x") + ggtitle("Density Distribution")
```

<img src="P8160_Project1_files/figure-gfm/unnamed-chunk-2-1.png" style="display: block; margin: auto;" />

## Generate Gompertz distribution using simsurv package

``` r
generate_gompertz = function(gamma, N, seed){
  set.seed(seed)
  #gen gompertz data
  covs <- data.frame(id = 1:N,
                    trt = stats::rbinom(N, 1, 0.5))
  dat <- simsurv(dist = "gompertz",
                 lambdas = 0.5, 
                 gammas = gamma, 
                 betas = c(trt = -0.5), 
                 x = covs, 
                 maxt = 5)
  dat <- merge(covs, dat)
  return(dat)
}
gompertz_dat <- generate_gompertz(1, 1000, 2023)

eventtime <- data.frame(gompertz_dat$eventtime)

# Visualization to validate the algorithm
eventtime %>% ggplot(aes(x = gompertz_dat.eventtime)) +
  geom_histogram(aes(y = after_stat(density)),
                 colour = 1, fill = "red", bins = 30) +
  geom_density() + xlab("x") + ggtitle("Density Distribution")
```

<img src="P8160_Project1_files/figure-gfm/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

## fit model

``` r
# fit Exponential
fit.exponential <- survreg(Surv(eventtime, status) ~ trt, data = gompertz_dat, dist = "exponential") 
summary(fit.exponential)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(eventtime, status) ~ trt, data = gompertz_dat, 
    ##     dist = "exponential")
    ##               Value Std. Error     z     p
    ## (Intercept) -0.0559     0.0440 -1.27 2e-01
    ## trt          0.2290     0.0633  3.62 3e-04
    ## 
    ## Scale fixed at 1 
    ## 
    ## Exponential distribution
    ## Loglik(model)= -1054.7   Loglik(intercept only)= -1061.2
    ##  Chisq= 13.1 on 1 degrees of freedom, p= 3e-04 
    ## Number of Newton-Raphson Iterations: 4 
    ## n= 1000

``` r
-fit.exponential$coefficients[-1]
```

    ##        trt 
    ## -0.2289529

``` r
# fit Weibull
fit.weibull <- survreg(Surv(eventtime, status) ~ trt, data = gompertz_dat, dist = "weibull") 
summary(fit.weibull)
```

    ## 
    ## Call:
    ## survreg(formula = Surv(eventtime, status) ~ trt, data = gompertz_dat, 
    ##     dist = "weibull")
    ##               Value Std. Error      z       p
    ## (Intercept)  0.0462     0.0295   1.57    0.12
    ## trt          0.2192     0.0414   5.30 1.2e-07
    ## Log(scale)  -0.4247     0.0260 -16.32 < 2e-16
    ## 
    ## Scale= 0.654 
    ## 
    ## Weibull distribution
    ## Loglik(model)= -943.7   Loglik(intercept only)= -957.5
    ##  Chisq= 27.62 on 1 degrees of freedom, p= 1.5e-07 
    ## Number of Newton-Raphson Iterations: 6 
    ## n= 1000

``` r
-fit.weibull$coefficients[-1] / fit.weibull$scale
```

    ##        trt 
    ## -0.3351914

``` r
#fit Cox model
fit.cox <- coxph(Surv(eventtime, status) ~ trt, data = gompertz_dat) 
summary(fit.cox)
```

    ## Call:
    ## coxph(formula = Surv(eventtime, status) ~ trt, data = gompertz_dat)
    ## 
    ##   n= 1000, number of events= 1000 
    ## 
    ##         coef exp(coef) se(coef)      z Pr(>|z|)    
    ## trt -0.40092   0.66970  0.06465 -6.201  5.6e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##     exp(coef) exp(-coef) lower .95 upper .95
    ## trt    0.6697      1.493      0.59    0.7602
    ## 
    ## Concordance= 0.547  (se = 0.009 )
    ## Likelihood ratio test= 38.46  on 1 df,   p=6e-10
    ## Wald test            = 38.46  on 1 df,   p=6e-10
    ## Score (logrank) test = 38.93  on 1 df,   p=4e-10

# Plot of Baseline Hazard Function

# Simulation for Gompertz

- N: sample size 100, 150, 200, 250, 300, 350, 400
- m: simulation time 1000
- $\beta$: true treatment effect to be -0.5
- $\lambda$: 0.5
- $\gamma$: 0.05, 1, 1.5

``` r
#write a fn to simulate gompertz data
sim_gompertz <- function(k, gamma=0.05, N){
  #generate gompertz data
  covs <- data.frame(id = 1:N,
                    trt = stats::rbinom(N, 1, 0.5))
  dat <- simsurv(dist = "gompertz",
                 lambdas = 0.5, 
                 gammas = gamma, 
                 betas = c(trt = -0.5), 
                 x = covs, 
                 maxt = 5)
  dat <- merge(covs, dat)
  #fit models
  fit.exponential <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "exponential")
  fit.weibull <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "weibull")
  fit.cox <- coxph(Surv(eventtime, status) ~ trt, data = dat)
  
  #extract beta
  result <- tibble(exp_beta = c(-fit.exponential$coefficients[-1]), 
                  weibull_beta = c(-fit.weibull$coefficients[-1])/fit.weibull$scale,
                  cox_beta = c(fit.cox$coefficients), 
                  dist = "gompertz",
                  beta = -0.5, 
                  gamma = gamma,
                  N = N)
  return(result)
  }

# Set seed for reproducibility
set.seed(2023)

# Create empty dataframe to store results
sim_gompertz_result1 <- data.frame()

# Simulate 1000 times
#gamma=0.05
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_gompertz(gamma = 0.05, N = n)
  sim_gompertz_result1 <- rbind(sim_gompertz_result1, sim_res)
}
}

gompertz_table1 <- sim_gompertz_result1 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
gompertz_table1
```

    ## # A tibble: 7 × 10
    ##       N mse_exp_ mse_weibull mse_cox var_exp var_weib…¹ var_cox bias_…² bias_w…³
    ##   <dbl>    <dbl>       <dbl>   <dbl>   <dbl>      <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100   0.0455      0.0508  0.0522  0.0454     0.0508  0.0522 0.0127  -0.00654
    ## 2   150   0.0280      0.0309  0.0317  0.0280     0.0308  0.0315 0.00617 -0.0104 
    ## 3   200   0.0212      0.0232  0.0239  0.0213     0.0230  0.0236 0.00259 -0.0141 
    ## 4   250   0.0174      0.0181  0.0186  0.0170     0.0181  0.0187 0.0222   0.00777
    ## 5   300   0.0142      0.0152  0.0157  0.0141     0.0152  0.0156 0.00958 -0.00574
    ## 6   350   0.0127      0.0133  0.0137  0.0124     0.0133  0.0137 0.0167   0.00243
    ## 7   400   0.0112      0.0119  0.0122  0.0109     0.0119  0.0122 0.0175   0.00264
    ## # … with 1 more variable: bias_cox <dbl>, and abbreviated variable names
    ## #   ¹​var_weibull, ²​bias_exp, ³​bias_weibull

``` r
#gamma=1
set.seed(2023)
sim_gompertz_result2 <- data.frame()
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_gompertz(gamma = 1, N = n)
  sim_gompertz_result2 <- rbind(sim_gompertz_result2, sim_res)
}
}

gompertz_table2 <- sim_gompertz_result2 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
gompertz_table2
```

    ## # A tibble: 7 × 10
    ##       N mse_e…¹ mse_w…² mse_cox var_exp var_w…³ var_cox bias_…⁴ bias_…⁵ bias_cox
    ##   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100  0.0655  0.0410  0.0494 0.0158  0.0359   0.0493   0.223  0.0712 -0.0129 
    ## 2   150  0.0568  0.0253  0.0289 0.00898 0.0207   0.0287   0.219  0.0683 -0.0144 
    ## 3   200  0.0535  0.0198  0.0220 0.00718 0.0157   0.0217   0.215  0.0640 -0.0197 
    ## 4   250  0.0573  0.0192  0.0177 0.00580 0.0124   0.0177   0.227  0.0824  0.00138
    ## 5   300  0.0535  0.0157  0.0147 0.00478 0.0105   0.0146   0.221  0.0723 -0.00981
    ## 6   350  0.0541  0.0152  0.0129 0.00415 0.00904  0.0129   0.223  0.0788 -0.00206
    ## 7   400  0.0536  0.0142  0.0110 0.00364 0.00819  0.0110   0.224  0.0779 -0.00323
    ## # … with abbreviated variable names ¹​mse_exp_, ²​mse_weibull, ³​var_weibull,
    ## #   ⁴​bias_exp, ⁵​bias_weibull

``` r
#gamma=1.5
set.seed(2023)
sim_gompertz_result3 <- data.frame()
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_gompertz(gamma = 1, N = n)
  sim_gompertz_result3  <- rbind(sim_gompertz_result3 , sim_res)
}
}

gompertz_table3 <- sim_gompertz_result3  %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
gompertz_table3
```

    ## # A tibble: 7 × 10
    ##       N mse_e…¹ mse_w…² mse_cox var_exp var_w…³ var_cox bias_…⁴ bias_…⁵ bias_cox
    ##   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100  0.0655  0.0410  0.0494 0.0158  0.0359   0.0493   0.223  0.0712 -0.0129 
    ## 2   150  0.0568  0.0253  0.0289 0.00898 0.0207   0.0287   0.219  0.0683 -0.0144 
    ## 3   200  0.0535  0.0198  0.0220 0.00718 0.0157   0.0217   0.215  0.0640 -0.0197 
    ## 4   250  0.0573  0.0192  0.0177 0.00580 0.0124   0.0177   0.227  0.0824  0.00138
    ## 5   300  0.0535  0.0157  0.0147 0.00478 0.0105   0.0146   0.221  0.0723 -0.00981
    ## 6   350  0.0541  0.0152  0.0129 0.00415 0.00904  0.0129   0.223  0.0788 -0.00206
    ## 7   400  0.0536  0.0142  0.0110 0.00364 0.00819  0.0110   0.224  0.0779 -0.00323
    ## # … with abbreviated variable names ¹​mse_exp_, ²​mse_weibull, ³​var_weibull,
    ## #   ⁴​bias_exp, ⁵​bias_weibull

# Simulation for Exponential

- N: sample size 100, 150, 200, 250, 300, 350, 400
- m: simulation time 1000
- $\beta$: true treatment effect to be -0.5
- $\lambda$: 0.5

``` r
#write a fn to simulate exponential data
sim_exp <- function(k, N){
  #generate exponential data
  covs <- data.frame(id = 1:N,
                    trt = stats::rbinom(N, 1, 0.5))
  dat <- simsurv(dist = "exponential",
                 lambdas = 0.5, 
                 betas = c(trt = -0.5), 
                 x = covs, 
                 maxt = 5)
  dat <- merge(covs, dat)
  #fit models
  fit.exponential <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "exponential")
  fit.weibull <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "weibull")
  fit.cox <- coxph(Surv(eventtime, status) ~ trt, data = dat)
  
  #extract beta
  result <- tibble(exp_beta = c(-fit.exponential$coefficients[-1]), 
                  weibull_beta = c(-fit.weibull$coefficients[-1])/fit.weibull$scale,
                  cox_beta = c(fit.cox$coefficients), 
                  dist = "exponential",
                  beta = -0.5, 
                  gamma = "default",
                  N = N)
  return(result)
  }

# Set seed for reproducibility
set.seed(2023)

# Create empty dataframe to store results
sim_exp_result1  <- data.frame()

# Simulate 1000 times
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_exp(N = n)
  sim_exp_result1 <- rbind(sim_exp_result1, sim_res)
}
}

exponential_table <- sim_exp_result1 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )

exponential_table
```

    ## # A tibble: 7 × 10
    ##       N mse_exp_ mse_weibull mse_cox var_exp var_wei…¹ var_cox bias_exp bias_w…²
    ##   <dbl>    <dbl>       <dbl>   <dbl>   <dbl>     <dbl>   <dbl>    <dbl>    <dbl>
    ## 1   100   0.0508      0.0539  0.0537  0.0508    0.0538  0.0537 -0.00585 -0.0121 
    ## 2   150   0.0315      0.0329  0.0326  0.0313    0.0326  0.0325 -0.0133  -0.0168 
    ## 3   200   0.0243      0.0251  0.0249  0.0241    0.0247  0.0247 -0.0162  -0.0196 
    ## 4   250   0.0194      0.0196  0.0196  0.0194    0.0196  0.0196  0.00304  0.00113
    ## 5   300   0.0157      0.0162  0.0162  0.0156    0.0160  0.0161 -0.00884 -0.0112 
    ## 6   350   0.0138      0.0141  0.0141  0.0138    0.0141  0.0141 -0.00330 -0.00478
    ## 7   400   0.0122      0.0126  0.0126  0.0122    0.0126  0.0126 -0.00188 -0.00388
    ## # … with 1 more variable: bias_cox <dbl>, and abbreviated variable names
    ## #   ¹​var_weibull, ²​bias_weibull

# Simulation for Weibull

- N: sample size 100, 150, 200, 250, 300, 350, 400
- m: simulation time 1000
- $\beta$: true treatment effect to be -0.5
- $\lambda$: 0.5
- $\gamma$: 0.05, 1, 1.5

``` r
#write a fn to simulate weibull data
sim_weibull <- function(k, gamma=0.05, N){
  #generate weibull data
  covs <- data.frame(id = 1:N,
                    trt = stats::rbinom(N, 1, 0.5))
  dat <- simsurv(dist = "weibull",
                 lambdas = 0.5, 
                 gammas = gamma, 
                 betas = c(trt = -0.5), 
                 x = covs, 
                 maxt = 5)
  dat <- merge(covs, dat)
  #fit models
  fit.exponential <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "exponential")
  fit.weibull <- survreg(Surv(eventtime, status) ~ trt, data = dat, dist = "weibull")
  fit.cox <- coxph(Surv(eventtime, status) ~ trt, data = dat)
  
  #extract beta
  result <- tibble(exp_beta = c(-fit.exponential$coefficients[-1]), 
                  weibull_beta = c(-fit.weibull$coefficients[-1])/fit.weibull$scale,
                  cox_beta = c(fit.cox$coefficients), 
                  dist = "weibull",
                  beta = -0.5, 
                  gamma = gamma,
                  N = N)
  return(result)
  }

# Set seed for reproducibility
set.seed(2023)

# Create empty dataframe to store results
sim_weibull_result1 <- data.frame()

# Simulate 1000 times
#gamma=0.05
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_weibull(gamma = 0.05, N = n)
  sim_weibull_result1 <- rbind(sim_weibull_result1, sim_res)
}
}

weibull_table1 <- sim_weibull_result1 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
weibull_table1
```

    ## # A tibble: 7 × 10
    ##       N  mse_exp_ mse_weibull mse_cox   var_exp var_w…¹ var_cox bias_…² bias_w…³
    ##   <dbl>     <dbl>       <dbl>   <dbl>     <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100 418805.        0.149   0.147  419199.    0.149   0.147   -4.98  -0.0228 
    ## 2   150   8027.        0.0805  0.0797   8020.    0.0798  0.0791  -3.91  -0.0279 
    ## 3   200     29.7       0.0673  0.0665     29.0   0.0664  0.0658  -0.812 -0.0304 
    ## 4   250     13.3       0.0499  0.0496     13.0   0.0499  0.0496  -0.596  0.00233
    ## 5   300      7.44      0.0417  0.0414      7.23  0.0417  0.0414  -0.467 -0.00956
    ## 6   350      1.45      0.0359  0.0356      1.30  0.0359  0.0356  -0.394 -0.00561
    ## 7   400      1.26      0.0315  0.0314      1.12  0.0313  0.0312  -0.380 -0.0138 
    ## # … with 1 more variable: bias_cox <dbl>, and abbreviated variable names
    ## #   ¹​var_weibull, ²​bias_exp, ³​bias_weibull

``` r
#gamma=1
sim_weibull_result2 <- data.frame()
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_weibull(gamma = 1, N = n)
  sim_weibull_result2 <- rbind(sim_weibull_result2, sim_res)
}
}

weibull_table2 <- sim_weibull_result1 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
weibull_table2
```

    ## # A tibble: 7 × 10
    ##       N  mse_exp_ mse_weibull mse_cox   var_exp var_w…¹ var_cox bias_…² bias_w…³
    ##   <dbl>     <dbl>       <dbl>   <dbl>     <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100 418805.        0.149   0.147  419199.    0.149   0.147   -4.98  -0.0228 
    ## 2   150   8027.        0.0805  0.0797   8020.    0.0798  0.0791  -3.91  -0.0279 
    ## 3   200     29.7       0.0673  0.0665     29.0   0.0664  0.0658  -0.812 -0.0304 
    ## 4   250     13.3       0.0499  0.0496     13.0   0.0499  0.0496  -0.596  0.00233
    ## 5   300      7.44      0.0417  0.0414      7.23  0.0417  0.0414  -0.467 -0.00956
    ## 6   350      1.45      0.0359  0.0356      1.30  0.0359  0.0356  -0.394 -0.00561
    ## 7   400      1.26      0.0315  0.0314      1.12  0.0313  0.0312  -0.380 -0.0138 
    ## # … with 1 more variable: bias_cox <dbl>, and abbreviated variable names
    ## #   ¹​var_weibull, ²​bias_exp, ³​bias_weibull

``` r
#gamma=1.5
sim_weibull_result3 <- data.frame()
for (n in c(100, 150, 200, 250, 300, 350, 400)) {
for (i in 1:1000) {
  sim_res <- sim_weibull(gamma = 1.5, N = n)
  sim_weibull_result3 <- rbind(sim_weibull_result3, sim_res)
}
}

weibull_table3 <- sim_weibull_result1 %>% 
  group_by(N) %>%
  summarize(mse_exp_ = mean((exp_beta+0.5)^2),
            mse_weibull = mean((weibull_beta+0.5)^2),
            mse_cox = mean((cox_beta+0.5)^2),
            var_exp = var(exp_beta),
            var_weibull = var(weibull_beta),
            var_cox = var(cox_beta),
            bias_exp = mean(exp_beta+0.5),
            bias_weibull = mean(weibull_beta+0.5),
            bias_cox = mean(cox_beta+0.5)
        
  )
weibull_table3
```

    ## # A tibble: 7 × 10
    ##       N  mse_exp_ mse_weibull mse_cox   var_exp var_w…¹ var_cox bias_…² bias_w…³
    ##   <dbl>     <dbl>       <dbl>   <dbl>     <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1   100 418805.        0.149   0.147  419199.    0.149   0.147   -4.98  -0.0228 
    ## 2   150   8027.        0.0805  0.0797   8020.    0.0798  0.0791  -3.91  -0.0279 
    ## 3   200     29.7       0.0673  0.0665     29.0   0.0664  0.0658  -0.812 -0.0304 
    ## 4   250     13.3       0.0499  0.0496     13.0   0.0499  0.0496  -0.596  0.00233
    ## 5   300      7.44      0.0417  0.0414      7.23  0.0417  0.0414  -0.467 -0.00956
    ## 6   350      1.45      0.0359  0.0356      1.30  0.0359  0.0356  -0.394 -0.00561
    ## 7   400      1.26      0.0315  0.0314      1.12  0.0313  0.0312  -0.380 -0.0138 
    ## # … with 1 more variable: bias_cox <dbl>, and abbreviated variable names
    ## #   ¹​var_weibull, ²​bias_exp, ³​bias_weibull

# Comparing the accuracy and efficiency of treatment effect ()

*k*: the number of independent data sets  

#### Bias

*Exponential Proportional-Hazards Model*:  
$$Bias=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{exp}^{(k)}}-\beta_{1})=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{exp}^{(k)}}-0.5)$$  

*Weibull Proportional-Hazards Model*:  
$$Bias=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{weibull}^{(k)}}-\beta_{1})=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{weibull}^{(k)}}-0.5)$$  

*Cox Proportional-Hazards Model*:  
$$Bias=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{cox}^{(k)}}-\beta_{1})=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{cox}^{(k)}}-0.5)$$  

#### Variance

*Exponential Proportional-Hazards Model*:  
$$Variance= \frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{exp}^{(k)}}-\beta_{1})^2=\frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{exp}^{(k)}}-0.5)^2$$  

*Weibull Proportional-Hazards Model*:  
$$Variance= \frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{weibull}^{(k)}}-\beta_{1})^2=\frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{weibull}^{(k)}}-0.5)^2$$  

*Cox Proportional-Hazards Model*:  
$$Variance= \frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{cox}^{(k)}}-\beta_{1})^2=\frac{1}{k-1}\sum_{i=1}^k\hat{(\beta_{cox}^{(k)}}-0.5)^2$$  

#### MSE

*Exponential Proportional-Hazards Model*:  
$$MSE= \frac{1}{k}\sum_{i=1}^k(\hat{\beta_{exp}^{(k)}}-\beta_{1})^2=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{exp}^{(k)}}-0.5)^2$$  

*Weibull Proportional-Hazards Model*:  
$$MSE= \frac{1}{k}\sum_{i=1}^k(\hat{\beta_{weibull}^{(k)}}-\beta_{1})^2=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{weibull}^{(k)}}-0.5)^2$$  

*Cox Proportional-Hazards Model*:  
$$MSE= \frac{1}{k}\sum_{i=1}^k(\hat{\beta_{cox}^{(k)}}-\beta_{1})^2=\frac{1}{k}\sum_{i=1}^k(\hat{\beta_{cox}^{(k)}}-0.5)^2$$  
