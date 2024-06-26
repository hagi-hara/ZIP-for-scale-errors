rm(list=ls())

library(cmdstanr)     # version 0.5.3
library(rstan)        # version 2.26.13
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.5.1
library(dplyr)        # version 1.1.1
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
library(tidyverse)    # version 2.0.0
library(bayesplot)    # version 1.10.0

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="classroom") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         session_z = scale(session, center=FALSE, scale=TRUE)[,1]-0.2276967,
         n_obj_c = n_obj-3)

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)
d %>% select(starts_with("session")) %>% distinct() %>% arrange(session)

# Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_sum,
              age = d$age_z,
              age2 = d$age_z^2,
              session = d$session_z,
              obj = d$n_obj_c,
              N_age_pred = length(age_pred$age_z),
              age_pred = age_pred$age_z,
              age2_pred = age_pred$age_z^2)


############################################################
##### ZIP Distribution with linear/quadratic functions #####
############################################################


################################################
# Model 1: bf(se_sum ~ age + age^2, 
#                 zi ~ age + age^2)
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: bf(se_sum ~ age + age^2, 
#                 zi ~ age        )
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model 3: bf(se_sum ~ age        , 
#                 zi ~ age + age^2)
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Model 4: bf(se_sum ~ age         , 
#                 zi ~ age         )
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model4.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4
print(fit4, pars = c("b_beta0", "b_beta1", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4, merge_chains = FALSE)
waic4 <- waic(log_lik)
waic4 %>% print(digits=2)


################################################
# Visualization
################################################
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "Analysis1", "Classroom")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit2_pred <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit2_pred %>% subset(category=="ber") %>% ggplot(aes(x=age))+
  geom_abline(intercept=0.5, slope=0, linetype="dashed", color="gray70", lwd=0.6)+
  geom_line(aes(y=MED, group=category), color="#1B9E77", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#1B9E77", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(13,41), breaks=seq(15,41,by=5))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Analysis1", "Classroom", "Figure_3a.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Count part
gp <- fit2_pred %>% subset(category=="poi") %>% ggplot(aes(x=age))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(13,41), breaks=seq(15,41,by=5))+
  scale_y_continuous(limits=c(0,12), breaks=0:12)+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Analysis1", "Classroom", "Figure_3b.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Overall ZIP
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 18, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit2_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit2_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(13,41), breaks=seq(15,41,by=5))+
  scale_y_continuous(limits=c(0,12), breaks=0:12)+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,5,10,15))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Analysis1", "Classroom", "Figure_3c.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit2_dist <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit2_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 4) %>% nrow() / nrow(d) # 95.0% of the actual data are within 95% of post predictive distribution [0, 4]

gp <- d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit2_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_classic()+
  labs(y="Density", x="The number of scale errors")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Analysis1", "Classroom", "Figure_fit2_PostPredictiveCheck.png"), plot=gp, dpi=350, width=3.0, height=4.4)  


################################################
# Age differences
################################################
fit2_agediff <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>%
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid") %>% 
  subset(age==26 | age ==27 | age==13 | age==41) %>% 
  select(-Parameter, -ageid, -age_z) %>% 
  pivot_wider(names_from = c(category, age), values_from = value) %>% 
  mutate(diff_Ber_26_minus_13 = ber_26 - ber_13,
         diff_Poi_26_miuns_13 = poi_26 - poi_13,
         diff_All_26_minus_13 = all_26 - all_13,
         diff_Ber_26_minus_41 = ber_26 - ber_41,
         diff_Poi_26_miuns_41 = poi_26 - poi_41,
         diff_All_26_minus_41 = all_26 - all_41,
         diff_All_27_minus_13 = all_27 - all_13,
         diff_All_27_minus_41 = all_27 - all_41,
         diff_Ber_41_minus_13 = ber_41 - ber_13) %>% 
  select(contains("diff")) %>% 
  pivot_longer(everything(), names_to = "diff", values_to = "value") %>% 
  group_by(diff) %>% 
  summarise(MED = median(value),
            lwr = quantile(value, prob=0.025),
            upr = quantile(value, prob=0.975))  
fit2_agediff %>% print(digits=2) 


############################################################
##### Poisson Distribution with linear/quadratic functions #
############################################################


################################################
# Model 5: se_sum ~ age + age^2
################################################

## For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)

## Make a data list for stan
dlist <- list(N = nrow(d),
              Y = d$se_sum,
              age = d$age_z,
              age2 = d$age_z^2,
              session = d$session_z,
              obj = d$n_obj_c,
              N_age_pred = length(age_pred$age_z),
              age_pred = age_pred$age_z,
              age2_pred = age_pred$age_z^2)

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model5.stan") %>%
  cmdstan_model() -> model5

fit5 <- model5$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model5")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit5$save_output_files(dir = output_dir, basename = "model5")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit5
print(fit5, pars = c("p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit5, merge_chains = FALSE)
waic5 <- waic(log_lik)
waic5 %>% print(digits=2)


################################################
# Model 6: se_sum ~ age
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis1", "Classroom", "model6.stan") %>%
  cmdstan_model() -> model6

fit6 <- model6$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "Analysis1", "Classroom", "model6")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit6$save_output_files(dir = output_dir, basename = "model6")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit6
print(fit6, pars = c("p_beta0", "p_beta1", "p_beta3", "p_beta4", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit6, merge_chains = FALSE)
waic6 <- waic(log_lik)
waic6 %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates,
      waic3$estimates,
      waic4$estimates,
      waic5$estimates,
      waic6$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:6), each = 3)) %>%
  cbind(Formula = rep(c("bf(se_sum ~ age + age^2, zi ~ age + age^2)",
                        "bf(se_sum ~ age + age^2, zi ~ age)",
                        "bf(se_sum ~ age, zi ~ age + age^2)",
                        "bf(se_sum ~ age, zi ~ age)",
                        "se_sum ~ age + age^2",
                        "se_sum ~ age"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic # the covariates are not written for simplicity

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)

# Correlation between theta and lambda in the selected model
ms <- rstan::extract(fit2)
N_mcmc <- length(ms$lp__)
r <- sapply(1:N_mcmc, function(i) cor(ms$lambda[i,], ms$theta[i,], method='spearman'))
quantile(r, prob=c(0.5, 0.025, 0.975))

# Convergence diagnostic check
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

ggs(fit2) %>% filter(str_detect(Parameter, pattern="b_")) %>% ggs_traceplot()
ggs(fit2) %>% filter(str_detect(Parameter, pattern="p_")) %>% ggs_traceplot()

ggs(fit2) %>% filter(str_detect(Parameter, pattern="b_")) %>% ggs_density()
ggs(fit2) %>% filter(str_detect(Parameter, pattern="p_")) %>% ggs_density()

neff_ratio(fit2)[1:8] 
mcmc_neff(neff_ratio(fit2)[1:8])
mcmc_neff_hist(neff_ratio(fit2)[1:8])+
  scale_x_continuous(limits=c(0,1))
