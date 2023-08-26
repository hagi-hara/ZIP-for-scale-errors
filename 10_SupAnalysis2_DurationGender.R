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
# The duration of scale errors: Inlab data
###################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & se_occ > 0 & se_dur!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="US" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-4,
         dur_c = dur-5)

d %>% group_by(groupid, country, session, n_obj, dur, n_obj_c, dur_c) %>% summarize(N=n())

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)

# Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_dur,
              age = d$age_z,
              age2 = d$age_z^2,
              session = d$session_num,
              obj = d$n_obj_c,
              gender = d$gender_num,
              N_age_pred = length(age_pred$age_z),
              age_pred = age_pred$age_z,
              age2_pred = age_pred$age_z^2)

################################################
# Model 1: se_dur ~ age + age^2
################################################

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "InLab", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "InLab", "model1")
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
here("Models", "SupAnalysis", "InLab", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "InLab", "model2")
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


################################################
# Visualization
################################################
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "SupAnalysis", "InLab")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit2_pred <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(ageid = as.numeric(str_sub(Parameter, start=8, end=-1))) %>% 
  left_join(d_age, by="ageid")

gp <- d %>% ggplot(aes(x=age))+
  geom_point(aes(y=se_dur), size=2, shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit2_pred, aes(y=MED), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit2_pred, aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,33), breaks=seq(15,33,by=3))+
  scale_y_continuous(limits=c(0,87), breaks=seq(0,80,by=10))+
  labs(y=expression(paste("Duration of scale errors (s)", sep="")), x="Age in  months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "InLab", "Figure_InLab_Dur.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit2_dist <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit2_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_dur >= 1.36365 & se_dur <= 125.017) %>% nrow() / nrow(d) # 98.9% of the actual data are within 95% of post predictive distribution

gp <- d %>% ggplot()+
  geom_histogram(aes(x=se_dur, y=..density..), bins=15, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit2_dist, aes(x=value, y=..density..), bins=15, color="white", fill="blue", alpha=0.3)+
  theme_classic()+
  labs(y="Density", x="Duration of scale errors")+
  scale_x_continuous(limits=c(0,200))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "InLab", "Figure_fit2_PostPredictiveCheck.png"), plot=gp, dpi=350, width=3.0, height=4.4)  


###################################################
# The duration of scale errors: Classroom data
###################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="classroom" & se_occ > 0 & se_dur!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         session_z = scale(session, center=FALSE, scale=TRUE)[,1]-0.184115,
         n_obj_c = n_obj-3)

d %>% group_by(groupid, n_obj, n_obj_c) %>% summarize(N=n())

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)
d %>% select(starts_with("session")) %>% distinct() %>% arrange(session)

# Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_dur,
              age = d$age_z,
              age2 = d$age_z^2,
              gender = d$gender_num,
              session = d$session_z,
              # obj = d$n_obj_c,
              N_age_pred = length(age_pred$age_z),
              age_pred = age_pred$age_z,
              age2_pred = age_pred$age_z^2)


################################################
# Model 1: se_dur ~ age + age^2
################################################

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Classroom", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Classroom", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("beta0", "beta1", "beta2", "beta3", "beta4", "sigma", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: se_dur ~ age
################################################

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Classroom", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Classroom", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("beta0", "beta1", "beta3", "beta4", "sigma", "lp__"), probs = c(0.5, 0.025, 0.975))

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

################################################
# Visualization
################################################
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "SupAnalysis", "Classroom")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit2_pred <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(ageid = as.numeric(str_sub(Parameter, start=8, end=-1))) %>% 
  left_join(d_age, by="ageid")

gp <- d %>% ggplot(aes(x=age))+
  geom_jitter(aes(y=se_dur), size=2, shape=16, color="gray60", alpha=0.4, height=0.1, width=0.1)+
  geom_line(data=fit2_pred, aes(y=MED), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit2_pred, aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(13.9,26.1), breaks=seq(14,26,by=3))+
  scale_y_continuous(limits=c(0,15.1), breaks=seq(0,15,by=3))+
  labs(y=expression(paste("Duration of scale errors (min)", sep="")), x="Age in  months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Classroom", "Figure_Classroom_Dur.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit2_dist <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit2_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_dur >= 0.28775498 & se_dur <= 9.04000100) %>% nrow() / nrow(d) # 98.1% of the actual data are within 95% of post predictive distribution

gp <- d %>% ggplot()+
  geom_histogram(aes(x=se_dur, y=..density..), bins=15, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit2_dist, aes(x=value, y=..density..), bins=15, color="white", fill="blue", alpha=0.3)+
  theme_classic()+
  labs(y="Density", x="Duration of scale errors")+
  scale_x_continuous(limits=c(0,40))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Classroom", "Figure_fit2_PostPredictiveCheck.png"), plot=gp, dpi=350, width=3.0, height=4.4)  
