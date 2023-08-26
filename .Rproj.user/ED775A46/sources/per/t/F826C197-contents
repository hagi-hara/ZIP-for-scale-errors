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
  subset(setting=="inlab") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="US" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000


####################################################################
##### ZIP Distribution with B-spline functions (Age, InLab) ########
####################################################################

################################################
# 5 base functions
################################################

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge). Change this value into 2 if you want make 6 base functions
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
              countUK = d$country_num_UK,
              countUS = d$country_num_US,
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit5$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1
print(fit1, pars = c("b_beta", "p_beta", 
                     "b_beta_countUK", "b_beta_countUS", "b_beta_session", "b_beta_obj", "b_beta_dur", "b_beta_gender",
                     "p_beta_countUK", "p_beta_countUS", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1, merge_chains = FALSE)
waic1 <- waic(log_lik)
waic1 %>% print(digits=2)


################################################
# Visualization
################################################
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age
output_fig_dir <- here("Figures", "SupAnalysis", "Bsplines")
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

# Logistic part
gp <- fit1_pred %>% subset(category=="ber") %>% ggplot(aes(x=age))+
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_InLab_Ber.png"), plot=gp, dpi=350, width=2.5, height=4.4)

# Count part
gp <- fit1_pred %>% subset(category=="poi") %>% ggplot(aes(x=age))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,6), breaks=0:6)+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_InLab_Poi.png"), plot=gp, dpi=350, width=2.5, height=4.4) 

# Overall ZIP
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 26, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit1_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit1_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,7), breaks=0:7)+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,5,10,20))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_InLab_All.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit5_dist <- ggs(fit1) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit1_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 5) %>% nrow() / nrow(d) # 98.2% of the actual data are within 95% of post predictive distribution

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit1_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
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


####################################################################
##### Poisson Distribution with B-spline functions (Age, InLab) ####
####################################################################

################################################
# 5 base functions
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
              countUK = d$country_num_UK,
              countUS = d$country_num_US,
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model2.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2
print(fit2, pars = c("p_beta", "p_beta_countUK", "p_beta_countUS", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2, merge_chains = FALSE)
waic2 <- waic(log_lik)
waic2 %>% print(digits=2)


################################################
# Visualization
################################################

age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "SupAnalysis", "Bsplines")
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

# Count part
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 26, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit2_pred %>% subset(category=="poi"), aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(data=fit2_pred %>% subset(category=="poi"), aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(15,37))+
  scale_y_continuous(limits=c(0,7), breaks=0:7)+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
  scale_size_continuous(breaks=c(1,5,10,20))+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=12, color="black", face="italic", hjust=0.2),
        legend.text=element_text(size=12, color="black"),
        legend.justification=c(1,0))
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_InLab_PoiOnly.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
# Post Predictive Check
fit2_dist <- ggs(fit2) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit2_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 4) %>% nrow() / nrow(d) # 95.4% of the actual data are within 95% of post predictive distribution

d %>% ggplot()+
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


########################################################################
##### ZIP Distribution with B-spline functions (Age, Classroom) ########
########################################################################

################################################
# 5 base functions
################################################

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge). Change this value into 2 if you want make 6 base functions
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
              session = d$session_z,
              obj = d$n_obj_c,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model3.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3
print(fit3, pars = c("b_beta", "p_beta", 
                     "b_beta_session", 
                     "p_beta_session", "p_beta_obj", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3, merge_chains = FALSE)
waic3 <- waic(log_lik)
waic3 %>% print(digits=2)


################################################
# Visualization
################################################
age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age
output_fig_dir <- here("Figures", "SupAnalysis", "Bsplines")
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_Cl_Ber.png"), plot=gp, dpi=350, width=2.5, height=4.4)

# Count part
gp <- fit3_pred %>% subset(category=="poi") %>% ggplot(aes(x=age))+
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_Cl_Poi.png"), plot=gp, dpi=350, width=2.5, height=4.4) 

# Overall ZIP
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 26, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit3_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit3_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_Cl_All.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit3_dist <- ggs(fit3) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit3_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 4) %>% nrow() / nrow(d) # 95.0% of the actual data are within 95% of post predictive distribution

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit3_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
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


########################################################################
##### Poisson Distribution with B-spline functions (Age, Classroom) ####
########################################################################

################################################
# 5 base functions
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
              session = d$session_z,
              obj = d$n_obj_c,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model4.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.9999, 
                      max_treedepth = 25,
                      seed = 456)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4
print(fit4, pars = c("p_beta", "p_beta_session", "p_beta_obj",  
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4, merge_chains = FALSE)
waic4 <- waic(log_lik)
waic4 %>% print(digits=2)


################################################
# Visualization
################################################

age_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age

output_fig_dir <- here("Figures", "SupAnalysis", "Bsplines")
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

# Count part
d %>% group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 18, MIN = 1

gp <- d %>%  group_by(age) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=age))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit4_pred %>% subset(category=="poi"), aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(data=fit4_pred %>% subset(category=="poi"), aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(13,41), breaks=seq(15,41,by=5))+
  scale_y_continuous(limits=c(0,12), breaks=0:12)+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Age in months")+
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_Cl_PoiOnly.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
# Post Predictive Check
fit4_dist <- ggs(fit4) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit4_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 4) %>% nrow() / nrow(d) # 95.0% of the actual data are within 95% of post predictive distribution

d %>% ggplot()+
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


########################################################################
##### ZIP Distribution with B-spline functions (Total Vocabulary) ######
########################################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & vcball!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         vcball_z = scale(vcball, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="US" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge)
x_min <- min(d$vcball_z)
x_max <- max(d$vcball_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$vcball_z, knots=knots, degree=3, intercept=TRUE)) # The number of basis functions = 1 (interior knots) + 4 (order of B-Spline) = 5

## Visualization of base functions
as.data.frame(t(bsMat)) %>% cbind(d$vcball_z) %>% 
  pivot_longer(!`d$vcball_z`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=`d$vcball_z`, y=value, color=basis_id))+
  geom_line()+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray60")+
  scale_x_continuous(limits=c(x_min, x_max))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()

## For prediction
vcb_pred <- d %>% select(starts_with("vcb")) %>% distinct() %>% arrange(vcball)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcball_z, knots=knots, degree=3, intercept=TRUE))


## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num, 
              N_age_pred = length(vcb_pred$vcball_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model5.stan") %>%
  cmdstan_model() -> model5

fit5 <- model5$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model5")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit5$save_output_files(dir = output_dir, basename = "model5")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit5
print(fit5, pars = c("b_beta", "p_beta",
                     "b_beta_countUK", "b_beta_session", "b_beta_obj", "b_beta_dur", "b_beta_gender",
                     "p_beta_countUK", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit5, merge_chains = FALSE)
waic5 <- waic(log_lik)
waic5 %>% print(digits=2)


################################################
# Visualization
################################################
vcb_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age
output_fig_dir <- here("Figures", "SupAnalysis", "Bsplines")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit5_pred <- ggs(fit5) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit5_pred %>% subset(category=="ber") %>% ggplot(aes(x=vcball))+
  geom_abline(intercept=0.5, slope=0, linetype="dashed", color="gray70", size=0.6)+
  geom_line(aes(y=MED, group=category), color="#1B9E77", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#1B9E77", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700), breaks=seq(0,700,by=150))+
  scale_y_continuous(limits=c(0,1))+
  labs(y=expression(paste("Mean proportion of producing scale errors (", italic(theta), ")", sep="")), x="Total vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_TVcb_Ber.png"), plot=gp, dpi=350, width=2.5, height=4.4)

# Count part
gp <- fit5_pred %>% subset(category=="poi") %>% ggplot(aes(x=vcball))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700), breaks=seq(0,700,by=150))+
  scale_y_continuous(limits=c(0,40))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Total vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_TVcb_Poi.png"), plot=gp, dpi=350, width=2.5, height=4.4) 

# Overall ZIP
d %>% group_by(vcball) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 4, MIN = 1

gp <- d %>%  group_by(vcball) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=vcball))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit5_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit5_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700), breaks=seq(0,700,by=150))+
  scale_y_continuous(limits=c(0,11), breaks=0:11)+
  labs(y=expression(paste("Mean number of scale errors (ZIP model)", sep="")), x="Total vocabulary size")+
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
ggsave(file = here("Figures", "SupAnalysis", "Bsplines", "Figure_BS_TVbc_ZIP.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Post Predictive Check
fit5_dist <- ggs(fit5) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

quantile(fit5_dist$value, c(0.0, 0.025, 0.5, 0.975, 1.0))
d %>% subset(se_sum >= 0 & se_sum <= 5) %>% nrow() / nrow(d) # 98.9% of the actual data are within 95% of post predictive distribution

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit5_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
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


########################################################################
##### ZIP Distribution with B-spline functions (Noun Vocabulary) #######
########################################################################

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

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge)
x_min <- min(d$vcb_z)
x_max <- max(d$vcb_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$vcb_z, knots=knots, degree=3, intercept=TRUE)) # The number of basis functions = 1 (interior knots) + 4 (order of B-Spline) = 5

## Visualization of base functions
as.data.frame(t(bsMat)) %>% cbind(d$vcb_z) %>% 
  pivot_longer(!`d$vcb_z`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=`d$vcb_z`, y=value, color=basis_id))+
  geom_line()+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray60")+
  scale_x_continuous(limits=c(x_min, x_max))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()

## For prediction
vcb_pred <- d %>% select(starts_with("vcb", ignore.case=FALSE)) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcb_z, knots=knots, degree=3, intercept=TRUE))

## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num, 
              N_age_pred = length(vcb_pred$vcb_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model5.stan") %>%
  cmdstan_model() -> model6

fit6 <- model6$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model6")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit6$save_output_files(dir = output_dir, basename = "model6")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit6
print(fit6, pars = c("b_beta", "p_beta",
                     "b_beta_countUK", "b_beta_session", "b_beta_obj", "b_beta_dur", "b_beta_gender",
                     "p_beta_countUK", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit6, merge_chains = FALSE)
waic6 <- waic(log_lik)
waic6 %>% print(digits=2)


########################################################################
##### ZIP Distribution with B-spline functions (Verb Vocabulary) #######
########################################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & noun!="NA") %>% 
  rename(Vcb = vcball,
         vcb = verb) %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         Vcb_z = scale(Vcb, center=TRUE, scale=TRUE)[,1],
         vcb_z = scale(vcb, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="USA" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge)
x_min <- min(d$vcb_z)
x_max <- max(d$vcb_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$vcb_z, knots=knots, degree=3, intercept=TRUE)) # The number of basis functions = 1 (interior knots) + 4 (order of B-Spline) = 5

## Visualization of base functions
as.data.frame(t(bsMat)) %>% cbind(d$vcb_z) %>% 
  pivot_longer(!`d$vcb_z`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=`d$vcb_z`, y=value, color=basis_id))+
  geom_line()+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray60")+
  scale_x_continuous(limits=c(x_min, x_max))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()

## For prediction
vcb_pred <- d %>% select(starts_with("vcb", ignore.case=FALSE)) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcb_z, knots=knots, degree=3, intercept=TRUE))


## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num, 
              N_age_pred = length(vcb_pred$vcb_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model5.stan") %>%
  cmdstan_model() -> model7

fit7 <- model7$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model7")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit7$save_output_files(dir = output_dir, basename = "model7")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit7
print(fit7, pars = c("b_beta", "p_beta",
                     "b_beta_countUK", "b_beta_session", "b_beta_obj", "b_beta_dur", "b_beta_gender",
                     "p_beta_countUK", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit7, merge_chains = FALSE)
waic7 <- waic(log_lik)
waic7 %>% print(digits=2)


########################################################################
##### ZIP Distribution with B-spline functions (Adj Vocabulary) ########
########################################################################

# Read the raw data
d <- read_csv("sedata.csv", col_types="ccccdcddccdddcddddddddddddddddddddddddddddddddd") %>% 
  subset(setting=="inlab" & noun!="NA") %>% 
  rename(Vcb = vcball,
         vcb = adj) %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         Vcb_z = scale(Vcb, center=TRUE, scale=TRUE)[,1],
         vcb_z = scale(vcb, center=TRUE, scale=TRUE)[,1],
         gender_num = case_when(gender=="f" ~ 1, gender=="m" ~ 0, TRUE ~ NA),
         country_num_UK = case_when(country=="UK" ~ 1, TRUE ~ 0),
         country_num_US = case_when(country=="USA" ~ 1, TRUE ~ 0),
         session_num = as.numeric(as.factor(session))-1,
         n_obj_c = n_obj-3,
         dur_c = dur-5)

## Make base functions
n_knots <- 1 # the number of interior knots (boundary knots will be added for each edge)
x_min <- min(d$vcb_z)
x_max <- max(d$vcb_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$vcb_z, knots=knots, degree=3, intercept=TRUE)) # The number of basis functions = 1 (interior knots) + 4 (order of B-Spline) = 5

## Visualization of base functions
as.data.frame(t(bsMat)) %>% cbind(d$vcb_z) %>% 
  pivot_longer(!`d$vcb_z`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=`d$vcb_z`, y=value, color=basis_id))+
  geom_line()+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray60")+
  scale_x_continuous(limits=c(x_min, x_max))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()

## For prediction
vcb_pred <- d %>% select(starts_with("vcb", ignore.case=FALSE)) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcb_z, knots=knots, degree=3, intercept=TRUE))


## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              countUK = d$country_num_UK, 
              session = d$session_num,
              obj = d$n_obj_c,
              dur = d$dur_c,
              gender = d$gender_num, 
              N_age_pred = length(vcb_pred$vcb_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "SupAnalysis", "Bsplines", "model5.stan") %>%
  cmdstan_model() -> model8

fit8 <- model8$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "SupAnalysis", "Bsplines", "model8")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit8$save_output_files(dir = output_dir, basename = "model8")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit8
print(fit8, pars = c("b_beta", "p_beta",
                     "b_beta_countUK", "b_beta_session", "b_beta_obj", "b_beta_dur", "b_beta_gender",
                     "p_beta_countUK", "p_beta_session", "p_beta_obj", "p_beta_dur", "p_beta_gender", 
                     "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit8, merge_chains = FALSE)
waic8 <- waic(log_lik)
waic8 %>% print(digits=2)

