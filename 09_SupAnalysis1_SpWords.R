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

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & noun_sp!="NA") %>% 
  mutate(gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

d %>% group_by(groupid, country, session, n_obj, dur) %>% summarize(N=n())

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000


################################################
# verb_sp = verb (only)
################################################

## Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_sum,
              noun_sp = d$noun_sp,
              verb_sp = d$verb_sp,   
              adj_sp = d$adj_sp,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "SpWords", "SpWords.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "SpWords", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "p_beta5",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# verb_sp = verb + preposition
################################################

## Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_sum,
              noun_sp = d$noun_sp,
              verb_sp = d$verb2_sp,   
              adj_sp = d$adj_sp,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "SpWords", "SpWords.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "SpWords",  "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "p_beta5",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

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
  cbind(Formula = rep(c("verb_sp = verb only",
                        "verb_sp = verb + preposition"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic # the covariates are not written for simplicity

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)
