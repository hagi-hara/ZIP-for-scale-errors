library(cmdstanr)     # version 0.3.0
library(rstan)        # version 2.21.2
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.4.1
library(dplyr)        # version 1.0.4
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
library(tidyverse)    # version 1.3.0

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcdccd") %>% subset(vcb!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         country_num = as.numeric(as.factor(country))-1,
         session_num = as.numeric(as.factor(session))-1,
         vcb_z = scale(vcb, center=TRUE, scale=TRUE)[,1])

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
age_pred <- d %>% select(starts_with("age")) %>% distinct() %>% arrange(age)

vcb_pred <- d %>% select(starts_with("vcb")) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=3) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]

# Make a data list for stan
dlist_age <- list(N = nrow(d),
                  Y = d$se_sum,
                  age = d$age_z,
                  age2 = d$age_z^2,
                  session = d$session_num,
                  N_age_pred = length(age_pred$age_z),
                  age_pred = age_pred$age_z,
                  age2_pred = age_pred$age_z^2)

dlist_vcb <- list(N = nrow(d),
                  Y = d$se_sum,
                  age = d$vcb_z,
                  age2 = d$vcb_z^2,
                  session = d$session_num,
                  N_age_pred = length(vcb_pred$vcb_z),
                  age_pred = vcb_pred$vcb_z,
                  age2_pred = vcb_pred$vcb_z^2)

################################################
# Model 1: bf(se_sum ~ age + age^2 + session, 
#                 zi ~ age + age^2 + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist_age,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: bf(se_sum ~ age + age^2 + session, 
#                 zi ~ age         + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist_age,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta3",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model 3: bf(se_sum ~ age         + session, 
#                 zi ~ age + age^2 + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist_age,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Model 4: bf(se_sum ~ age         + session, 
#                 zi ~ age         + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model4.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist_age,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4
print(fit4, pars = c("b_beta0", "b_beta1", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4, merge_chains = FALSE)
waic4 <- waic(log_lik)
waic4 %>% print(digits=2)


################################################
# Model 5: bf(se_sum ~ vcb + vcb^2 + session, 
#                 zi ~ vcb + vcb^2 + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model1.stan") %>%
  cmdstan_model() -> model5

fit5 <- model5$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model5")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit5$save_output_files(dir = output_dir, basename = "model5")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit5
print(fit5, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit5, merge_chains = FALSE)
waic5 <- waic(log_lik)
waic5 %>% print(digits=2)


################################################
# Model 6: bf(se_sum ~ vcb + vcb^2 + session, 
#                 zi ~ vcb         + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model2.stan") %>%
  cmdstan_model() -> model6

fit6 <- model6$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model6")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit6$save_output_files(dir = output_dir, basename = "model6")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit6
print(fit6, pars = c("b_beta0", "b_beta1", "b_beta3",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit6, merge_chains = FALSE)
waic6 <- waic(log_lik)
waic6 %>% print(digits=2)


################################################
# Model 7: bf(se_sum ~ vcb         + session, 
#                 zi ~ vcb + vcb^2 + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model3.stan") %>%
  cmdstan_model() -> model7

fit7 <- model7$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model7")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit7$save_output_files(dir = output_dir, basename = "model7")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit7
print(fit7, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit7, merge_chains = FALSE)
waic7 <- waic(log_lik)
waic7 %>% print(digits=2)


################################################
# Model 8: bf(se_sum ~ vcb         + session, 
#                 zi ~ vcb         + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip_vcb", "model4.stan") %>%
  cmdstan_model() -> model8

fit8 <- model8$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_vcb", "model8")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit8$save_output_files(dir = output_dir, basename = "model8")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit8
print(fit8, pars = c("b_beta0", "b_beta1", "b_beta3", 
                     "p_beta0", "p_beta1", "p_beta3", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit8, merge_chains = FALSE)
waic8 <- waic(log_lik)
waic8 %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates,
      waic3$estimates,
      waic4$estimates,
      waic5$estimates,
      waic6$estimates,
      waic7$estimates,
      waic8$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:8), each = 3)) %>%
  cbind(Formula = rep(c("bf(se_sum ~ age + age^2, zi ~ age + age^2)",
                        "bf(se_sum ~ age + age^2, zi ~ age)",
                        "bf(se_sum ~ age, zi ~ age + age^2)",
                        "bf(se_sum ~ age, zi ~ age)",
                        "bf(se_sum ~ vcb + vcb^2, zi ~ vcb + vcb^2)",
                        "bf(se_sum ~ vcb + vcb^2, zi ~ vcb)",
                        "bf(se_sum ~ vcb, zi ~ vcb + vcb^2)",
                        "bf(se_sum ~ vcb, zi ~ vcb)"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic # the covariates are not written for simplicity

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)

# Correlation between theta and lambda in the selected model
ms <- rstan::extract(fit3)
N_mcmc <- length(ms$lp__)
r <- sapply(1:N_mcmc, function(i) cor(ms$lambda[i,], ms$theta[i,], method='spearman'))
quantile(r, prob=c(0.5, 0.025, 0.975))

ms <- rstan::extract(fit7)
N_mcmc <- length(ms$lp__)
r <- sapply(1:N_mcmc, function(i) cor(ms$lambda[i,], ms$theta[i,], method='spearman'))
quantile(r, prob=c(0.5, 0.025, 0.975))


################################################
# Visualization
################################################

# Model 3
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit3_pred <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit3_pred %>% subset(category=="ber") %>% ggplot(aes(x=age))+
  geom_abline(intercept=0.5, slope=0, linetype="dashed", color="gray70", size=0.6)+
  geom_line(aes(y=MED, group=category), color="#1B9E77", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#1B9E77", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,35))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_4_1.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Count part
gp <- fit3_pred %>% subset(category=="poi") %>% ggplot(aes(x=age))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,35))+
  scale_y_continuous(limits=c(0,6.1))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_4_2.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Overall ZIP
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 11, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit3_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit3_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,35))+
  scale_y_continuous(limits=c(0,7))+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,4,7,11))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Figure_4_3.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit3_dist <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit3_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()

# Model 7
vcb_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit7_pred <- ggs(fit7) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit7_pred %>% subset(category=="ber") %>% ggplot(aes(x=vcb))+
  geom_abline(intercept=0.5, slope=0, linetype="dashed", color="gray70", size=0.6)+
  geom_line(aes(y=MED, group=category), color="#1B9E77", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#1B9E77", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Total vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_4_4.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Count part
gp <- fit7_pred %>% subset(category=="poi") %>% ggplot(aes(x=vcb))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700))+
  scale_y_continuous(limits=c(0,6.1))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Total vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_4_5.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Overall ZIP
d %>% group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 3, MIN = 1

gp <- d %>%  group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=vcb))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit7_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit7_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700))+
  scale_y_continuous(limits=c(0,7))+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Total vocabulary size")+
  scale_size_continuous(breaks=c(1,2,3))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Figure_4_6.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit7_dist <- ggs(fit7) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit7_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()


