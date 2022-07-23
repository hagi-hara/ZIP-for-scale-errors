library(cmdstanr)     # version 0.3.0
library(rstan)        # version 2.21.2
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.4.1
library(dplyr)        # version 1.0.4
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
library(splines2)     # version 0.4.5
library(tidyverse)    # version 1.3.0

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcdccd") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         country_num = as.numeric(as.factor(country))-1,
         session_num = as.numeric(as.factor(session))-1)

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)

# Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_sum,
              age = d$age_z,
              age2 = d$age_z^2,
              country = d$country_num,
              session = d$session_num,
              N_age_pred = length(age_pred$age_z),
              age_pred = age_pred$age_z,
              age2_pred = age_pred$age_z^2)


################################################
# Model 1: se_sum ~ age + age^2 + country + session
################################################

## MCMC sampling with cmdstanr
here("Models", "Poisson", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Poisson", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: se_sum ~ age + country + session
################################################

## MCMC sampling with cmdstanr
here("Models", "Poisson", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Poisson", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("p_beta0", "p_beta1", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model 3. B-splines (5 base functions)
################################################

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge)
x_min <- min(d$age_z)
x_max <- max(d$age_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$age_z, knots=knots, degree=3, intercept=TRUE)) # The number of basis functions = 1 (interior knots) + 4 (order of B-Spline) = 5

## Visualization of base functions
as.data.frame(t(bsMat)) %>% cbind(d$age_z) %>% 
  pivot_longer(!`d$age_z`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=`d$age_z`, y=value, color=basis_id))+
  geom_line()+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray60")+
  scale_x_continuous(limits=c(x_min, x_max))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()

## For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)
bsMat_pred <- t(bSpline(age_pred$age_z, knots=knots, degree=3, intercept=TRUE))

## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              country = d$country_num,
              session = d$session_num,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "Poisson", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Poisson", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("p_beta", "p_beta_country", "p_beta_session", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates,
      waic3$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:3), each = 3)) %>%
  cbind(Formula = rep(c("se_sum ~ age + age^2",
                        "se_sum ~ age",
                        "B-splines (5 base functions)"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic # the covariates are not written for simplicity

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)


################################################
# Visualization
################################################

# Model 1
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit1_pred <- ggs(fit1) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Count part
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 22, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit1_pred %>% subset(category=="poi"), aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(data=fit1_pred %>% subset(category=="poi"), aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,7))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,8,15,22))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Figure_S3_1.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit1_dist <- ggs(fit1) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit1_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()


# Model 3
fit3_pred <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Count part
gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit3_pred %>% subset(category=="poi"), aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(data=fit3_pred %>% subset(category=="poi"), aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,7))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,8,15,22))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Figure_S3_2.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit3_dist <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit3_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()

