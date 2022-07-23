library(cmdstanr)     # version 0.3.0
library(rstan)        # version 2.21.2
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.4.1
library(dplyr)        # version 1.0.4
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
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
# Model 1: bf(se_sum ~ age + age^2 + country + session, 
#                 zi ~ age + age^2 + country + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", "b_beta4",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: bf(se_sum ~ age + age^2 + country + session, 
#                 zi ~ age         + country + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta3", "b_beta4",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model 3: bf(se_sum ~ age         + country + session, 
#                 zi ~ age + age^2 + country + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", "b_beta4",
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Model 4: bf(se_sum ~ age         + country + session, 
#                 zi ~ age         + country + session)
################################################

## MCMC sampling with cmdstanr
here("Models", "Zip", "model4.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4
print(fit4, pars = c("b_beta0", "b_beta1", "b_beta3", "b_beta4",
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4, merge_chains = FALSE)
waic4 <- waic(log_lik)
waic4 %>% print(digits=2)



################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates,
      waic3$estimates,
      waic4$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:4), each = 3)) %>%
  cbind(Formula = rep(c("bf(se_sum ~ age + age^2, zi ~ age + age^2)",
                        "bf(se_sum ~ age + age^2, zi ~ age)",
                        "bf(se_sum ~ age, zi ~ age + age^2)",
                        "bf(se_sum ~ age, zi ~ age)"), each = 3)) %>%
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


################################################
# Visualization
################################################
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
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_3_1.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Count part
gp <- fit3_pred %>% subset(category=="poi") %>% ggplot(aes(x=age))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,4.2))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_3_2.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Overall ZIP
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 22, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit3_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit3_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,7))+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Age in months")+
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
ggsave(file = here("Figures", "Figure_3_3.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit3_dist <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit3_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()


################################################
# Age differences
################################################
fit3_agediff <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>%
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid") %>% 
  subset(age==18 | age==22) %>% 
  select(-Parameter, -ageid, -age_z) %>% 
  pivot_wider(names_from = c(category, age), values_from = value) %>% 
  mutate(diff_Ber_22_minus_18 = ber_22 - ber_18,
         diff_Poi_22_miuns_18 = poi_22 - poi_18,
         diff_All_22_minus_18 = all_22 - all_18) %>% 
  select(contains("diff")) %>% 
  pivot_longer(everything(), names_to = "diff", values_to = "value") %>% 
  group_by(diff) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025))  
fit3_agediff %>% print(digits=2) 


################################################
# ZIP modeling using brms package (Example)
################################################
library(brms)
formula <- bf(se_sum ~ age_z + I(age_z^2), zi ~ age_z + I(age_z^2))
get_prior(formula,
          data = d,
          family=zero_inflated_poisson())
## MCMC sampling with brms
fit1_brms <- brm(formula,
                 data = d,
                 family=zero_inflated_poisson(),
                 seed = 321,
                 iter = iter_set,
                 warmup = warmup_set,
                 chains = 4,
                 control = list(max_treedepth = 15, adapt_delta = 0.995),
                 prior=c(set_prior("student_t(3, 0, 1)", class="b")),
                 backend="cmdstanr")

## Confirm priors and convergence
prior_summary(fit1_brms)
plot(fit1_brms, combo=c("trace", "dens_overlay"))

## Results
fit1_brms %>% summary()
brms::waic(fit1_brms) %>% print(digits=2)


################################################
# Poisson modeling using brms package (For reference)
################################################
## MCMC sampling with brms
fit2_brms <- brm(se_sum ~ age_z + I(age_z^2),
                 data = d,
                 family=poisson(),
                 seed = 321,
                 iter = iter_set,
                 warmup = warmup_set,
                 chains = 4,
                 control = list(max_treedepth = 15, adapt_delta = 0.995),
                 prior=c(set_prior("", class="b")),
                 backend="cmdstanr")

## Confirm priors and convergence
prior_summary(fit2_brms)
plot(fit2_brms, combo=c("trace", "dens_overlay"))

## Results
fit2_brms %>% summary()
brms::waic(fit2_brms) %>% print(digits=2)
