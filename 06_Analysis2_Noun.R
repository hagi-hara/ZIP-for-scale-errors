rm(list=ls())

library(cmdstanr)     # version 0.5.3
library(rstan)        # version 2.26.13
options(mc.cores = parallel::detectCores())
library(loo)          # version 2.5.1
library(dplyr)        # version 1.1.1
library(here)         # version 1.0.1
library(ggmcmc)       # version 1.5.1.1
library(tidyverse)    # version 2.0.0

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & noun!="NA") %>% 
  rename(Vcb = vcball,
         vcb = noun) %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         Vcb_z = scale(Vcb, center=TRUE, scale=TRUE)[,1],
         vcb_z = scale(vcb, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="USA" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

d %>% group_by(groupid, country, session, n_obj, dur) %>% summarize(N=n())

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

# For prediction
vcb_pred <- d %>% select(starts_with("vcb", ignore.case=FALSE)) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=3) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]

# Make a data list for stan
dlist_vcb <- list(N = nrow(d),
                  Y = d$se_sum,
                  age = d$vcb_z,
                  age2 = d$vcb_z^2,
                  countUK = d$country_num_UK, 
                  session = d$session_num,
                  obj = d$n_obj_c,
                  dur = d$dur_c,
                  gender = d$gender_num,                 
                  N_age_pred = length(vcb_pred$vcb_z),
                  age_pred = vcb_pred$vcb_z,
                  age2_pred = vcb_pred$vcb_z^2)


################################################
# Model 1: bf(se_sum ~ vcb + vcb^2, 
#                 zi ~ vcb + vcb^2)
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis2", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Analysis2", "Noun", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", "b_beta4", "b_beta5", "b_beta6", "b_beta7",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "p_beta5", "p_beta6", "p_beta7",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Model 2: bf(se_sum ~ vcb + vcb^2, 
#                 zi ~ vcb        )
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis2", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Analysis2", "Noun", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("b_beta0", "b_beta1", "b_beta3", "b_beta4", "b_beta5", "b_beta6", "b_beta7",
                     "p_beta0", "p_beta1", "p_beta2", "p_beta3", "p_beta4", "p_beta5", "p_beta6", "p_beta7",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Model 3: bf(se_sum ~ vcb        , 
#                 zi ~ vcb + vcb^2)
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis2", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Analysis2", "Noun", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", "b_beta4", "b_beta5", "b_beta6", "b_beta7",
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", "p_beta5", "p_beta6", "p_beta7",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Model 4: bf(se_sum ~ vcb        , 
#                 zi ~ vcb        )
################################################

## MCMC sampling with cmdstanr
here("Models", "Analysis2", "model4.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist_vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Analysis2", "Noun", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4
print(fit4, pars = c("b_beta0", "b_beta1", "b_beta3", "b_beta4", "b_beta5", "b_beta6", "b_beta7",
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", "p_beta5", "p_beta6", "p_beta7",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4, merge_chains = FALSE)
waic4 <- waic(log_lik)
waic4 %>% print(digits=2)


################################################
# Visualization
################################################

# Model 4
vcb_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "Analysis2", "Noun")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit4_pred <- ggs(fit4) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit4_pred %>% subset(category=="ber") %>% ggplot(aes(x=vcb))+
  geom_abline(intercept=0.5, slope=0, linetype="dashed", color="gray70", size=0.6)+
  geom_line(aes(y=MED, group=category), color="#1B9E77", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#1B9E77", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,260), breaks=seq(0,260,by=50))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Noun vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Analysis2", "Noun", "Figure_Noun_Ber.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Count part
gp <- fit4_pred %>% subset(category=="poi") %>% ggplot(aes(x=vcb))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,260), breaks=seq(0,260,by=50))+
  scale_y_continuous(limits=c(0,9.5), breaks=0:9)+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Noun vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Analysis2", "Noun", "Figure_Noun_Poi.png"), plot=gp, dpi=350, width=2.5, height=4.4)  

# Overall ZIP
d %>% group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 3, MIN = 1

gp <- d %>%  group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=vcb))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit4_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit4_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,260), breaks=seq(0,260,by=50))+
  scale_y_continuous(limits=c(0,9.5), breaks=0:9)+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Noun vocabulary size")+
  scale_size_continuous(breaks=c(1,2,3,4))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "Analysis2", "Noun", "Figure_Noun_ZIP.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit4_dist <- ggs(fit4) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit4_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 4) %>% nrow() / nrow(d) # 97.5% of the actual data are within 95% of post predictive distribution, [0,4]

gp <- d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit4_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
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
ggsave(file = here("Figures", "Analysis2", "Noun", "Figure_fit3_PostPredictiveCheck.png"), plot=gp, dpi=350, width=3.0, height=4.4)  


################################################
# Age differences
################################################
fit4_agediff <- ggs(fit4) %>%
  filter(str_detect(Parameter, pattern="_pred")) %>%
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>%
  left_join(d_age, by="ageid") %>%
  subset(vcb==0 | vcb==258) %>%
  select(-Parameter, -ageid, -vcb_z) %>%
  pivot_wider(names_from = c(category, vcb), values_from = value) %>%
  mutate(diff_Ber_258_minus_0 = ber_258 - ber_0,
         diff_Poi_258_miuns_0 = poi_258 - poi_0,
         diff_All_258_minus_0 = all_258 - all_0) %>%
  select(contains("diff")) %>%
  pivot_longer(everything(), names_to = "diff", values_to = "value") %>%
  group_by(diff) %>%
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025))
fit4_agediff %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1$estimates,
      waic2$estimates,
      waic3$estimates,
      waic4$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:4), each = 3)) %>%
  cbind(Formula = rep(c("bf(se_sum ~ vcb + vcb^2, zi ~ vcb + vcb^2)",
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
ms <- rstan::extract(fit4)
N_mcmc <- length(ms$lp__)
r <- sapply(1:N_mcmc, function(i) cor(ms$lambda[i,], ms$theta[i,], method='spearman'))
quantile(r, prob=c(0.5, 0.025, 0.975))


################################################
# Model 3:REFERENCE (Total Vocabulary Size) 
#         bf(se_sum ~ vcb        , 
#                 zi ~ vcb + vcb^2)
################################################

## For prediction
Vcb_pred <- d %>% select(starts_with("Vcb", ignore.case=FALSE)) %>% distinct() %>% arrange(Vcb)
seq <- seq(1, nrow(Vcb_pred), by=3) %>% append(nrow(Vcb_pred))
Vcb_pred <- Vcb_pred[seq,]

## Make a data list for stan
dlist_Vcb <- list(N = nrow(d),
                  Y = d$se_sum,
                  age = d$Vcb_z,
                  age2 = d$Vcb_z^2,
                  countUK = d$country_num_UK, 
                  session = d$session_num,
                  obj = d$n_obj_c,
                  dur = d$dur_c,
                  gender = d$gender_num,                 
                  N_age_pred = length(Vcb_pred$Vcb_z),
                  age_pred = Vcb_pred$Vcb_z,
                  age2_pred = Vcb_pred$Vcb_z^2)

## MCMC sampling with cmdstanr
here("Models", "Analysis2", "model3.stan") %>%
  cmdstan_model() -> model5

fit5 <- model5$sample(data = dlist_Vcb,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Analysis2", "TotalVcb2", "model5")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit5$save_output_files(dir = output_dir, basename = "model5")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit5
print(fit5, pars = c("b_beta0", "b_beta1", "b_beta2", "b_beta3", "b_beta4", "b_beta5", "b_beta6", "b_beta7",
                     "p_beta0", "p_beta1", "p_beta3", "p_beta4", "p_beta5", "p_beta6", "p_beta7",
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit5, merge_chains = FALSE)
waic5 <- waic(log_lik)
waic5 %>% print(digits=2)
