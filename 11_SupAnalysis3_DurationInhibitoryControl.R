rm(list=ls())

library(cmdstanr)     # version 0.5.3
library(rstan)        # version 2.26.13
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.5.1
library(dplyr)        # version 1.1.1
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
library(tidyverse)    # version 2.0.0
library(splines2)     # version 0.5.0

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000


###################################################
# The Effect of Inhibitory control on Scale Errors
###################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & se_occ > 0 & ec!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         ec_z = scale(ec, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

d %>% group_by(groupid, country, session, n_obj, dur) %>% summarize(N=n())

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)
ec_pred <- d %>% select(starts_with("ec", ignore.case=FALSE)) %>% distinct() %>% arrange(ec)
seq <- seq(1, nrow(ec_pred), by=3) %>% append(nrow(ec_pred))
ec_pred <- ec_pred[seq,]


# Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_dur,
              age = d$age_z,
              age2 = d$age_z^2,
              ec = d$ec_z,
              session = d$session_num,
              gender = d$gender_num)


################################################
# Model 1: se_dur ~ age + age^2
################################################

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "InhibitoryControl", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "InhibitoryControl", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "sigma", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: se_dur ~ age
################################################

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "InhibitoryControl", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "InhibitoryControl", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("beta0", "beta1", "beta3", "beta4", "beta5", "sigma", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:2), each = 3)) %>%
  cbind(Formula = rep(c("se_dur ~ age + age^2",
                        "se_dur ~ age"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic # the covariates are not written for simplicity

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)
