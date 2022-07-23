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
d <- read_csv("sedata.csv", col_types="ccccdcdccd") %>% subset(vcb!="NA") %>% 
  mutate(age_z = scale(age, center=TRUE, scale=TRUE)[,1],
         country_num = as.numeric(as.factor(country))-1,
         session_num = as.numeric(as.factor(session))-1,
         vcb_z = scale(vcb, center=TRUE, scale=TRUE)[,1])

# MCMC sampling setting
iter_set <- 15000
warmup_set <- 2000

################################################
# Model 1. 5 base functions (Age)
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
              session = d$session_num,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "Zip_bspline_vcb", "model1.stan") %>%
  cmdstan_model() -> model1

fit1 <- model1$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_bspline_vcb", "model1")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit1$save_output_files(dir = output_dir, basename = "model1")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit1.bs
print(fit1.bs, pars = c("b_beta", "p_beta", "b_beta_session", "p_beta_session", 
                        "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit1.bs, merge_chains = FALSE)
waic1.bs <- waic(log_lik)
waic1.bs %>% print(digits=2)


################################################
# Model 2. 6 base functions (Age)
################################################

## Make base functions
n_knots <- 2 
x_min <- min(d$age_z)
x_max <- max(d$age_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$age_z, knots=knots, degree=3, intercept=TRUE))

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
              session = d$session_num,
              N_age_pred = length(age_pred$age_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "Zip_bspline_vcb", "model1.stan") %>%
  cmdstan_model() -> model2

fit2 <- model2$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_bspline_vcb", "model2")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit2$save_output_files(dir = output_dir, basename = "model2")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit2.bs
print(fit2.bs, pars = c("b_beta", "p_beta", "b_beta_session", "p_beta_session", 
                        "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit2.bs, merge_chains = FALSE)
waic2.bs <- waic(log_lik)
waic2.bs %>% print(digits=2)


################################################
# Model 3. 5 base functions (Vocabulary)
################################################

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
vcb_pred <- d %>% select(starts_with("vcb")) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcb_z, knots=knots, degree=3, intercept=TRUE))


## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              session = d$session_num,
              N_age_pred = length(vcb_pred$vcb_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "Zip_bspline_vcb", "model1.stan") %>%
  cmdstan_model() -> model3

fit3 <- model3$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_bspline_vcb", "model3")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit3$save_output_files(dir = output_dir, basename = "model3")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit3.bs
print(fit3.bs, pars = c("b_beta", "p_beta", "b_beta_session", "p_beta_session", 
                        "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit3.bs, merge_chains = FALSE)
waic3.bs <- waic(log_lik)
waic3.bs %>% print(digits=2)


################################################
# Model 4. 6 base functions (Vocabulary)
################################################

## Make base functions
n_knots <- 2 
x_min <- min(d$vcb_z)
x_max <- max(d$vcb_z)
# knots <- quantile(d$age_z, seq(0,1,length=n_knots+2)) %>% as.vector() # interior knots + boundary knots of (mix_n, max_x)
knots <- seq(x_min, x_max, length=n_knots+2)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(d$vcb_z, knots=knots, degree=3, intercept=TRUE))

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
vcb_pred <- d %>% select(starts_with("vcb")) %>% distinct() %>% arrange(vcb)
seq <- seq(1, nrow(vcb_pred), by=10) %>% append(nrow(vcb_pred))
vcb_pred <- vcb_pred[seq,]
bsMat_pred <- t(bSpline(vcb_pred$vcb_z, knots=knots, degree=3, intercept=TRUE))

## Make a data list for stan
dlist <- list(N_row = nrow(d),
              N_basis = nrow(bsMat),
              Y = d$se_sum,
              B = bsMat,
              session = d$session_num,
              N_age_pred = length(vcb_pred$vcb_z),
              B_pred = bsMat_pred)

## MCMC sampling with cmdstanr
here("Models", "Zip_bspline_vcb", "model1.stan") %>%
  cmdstan_model() -> model4

fit4 <- model4$sample(data = dlist,
                      chains = 4,
                      iter_warmup = warmup_set,
                      iter_sampling = iter_set,
                      adapt_delta = 0.995, 
                      max_treedepth = 15,
                      seed = 321)

output_dir <- here("MCMCSamples", "Zip_bspline_vcb", "model4")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}
fit4$save_output_files(dir = output_dir, basename = "model4")

## get MCMC samples with stanfit object format
str_c(output_dir,
      list.files(output_dir, pattern = ".csv"), sep = "/") %>% 
  read_stan_csv() -> fit4.bs
print(fit4.bs, pars = c("b_beta", "p_beta", "b_beta_session", "p_beta_session", 
                        "lp__"), probs = c(0.5, 0.025, 0.975))

## calculate WAIC
log_lik <- extract_log_lik(fit4.bs, merge_chains = FALSE)
waic4.bs <- waic(log_lik)
waic4.bs %>% print(digits=2)


################################################
# Model Selection
################################################
rbind(waic1.bs$estimates,
      waic2.bs$estimates,
      waic3.bs$estimates,
      waic4.bs$estimates) %>% 
  cbind(Model = rep(str_c("model", 1:4), each = 3)) %>%
  cbind(Formula = rep(c("base functions = 5 (Age)",
                        "base functions = 6 (Age)",
                        "base functions = 5 (Vocabulary)",
                        "base functions = 6 (Vocabulary)"), each = 3)) %>%
  cbind(Measure = c("elpd_waic", "p_waic", "waic")) %>% 
  as_tibble() -> df_waic

df_waic %>% 
  filter(Measure == "waic") %>%
  relocate(Model, Formula, Measure) %>% 
  arrange(Estimate)

# Correlation between theta and lambda in the selected model
ms <- rstan::extract(fit3.bs)
N_mcmc <- length(ms$lp__)
r <- sapply(1:N_mcmc, function(i) cor(ms$lambda[i,], ms$theta[i,], method='spearman'))
quantile(r, prob=c(0.5, 0.025, 0.975))


################################################
# Visualization
################################################
vcb_pred %>% tibble %>% mutate(ageid = row_number()) -> d_age
output_fig_dir <- here("Figures")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

fit3_pred <- ggs(fit3.bs) %>% 
  filter(str_detect(Parameter, pattern="_pred")) %>% 
  group_by(Parameter) %>% 
  summarise(MED = median(value),
            upr = quantile(value, prob=0.975),
            lwr = quantile(value, prob=0.025)) %>% 
  mutate(category = str_sub(Parameter, start=8, end=10),
         ageid = as.numeric(str_sub(Parameter, start=12, end=-1))) %>% 
  left_join(d_age, by="ageid")

# Logistic part
gp <- fit3_pred %>% subset(category=="ber") %>% ggplot(aes(x=vcb))+
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
ggsave(file = here("Figures", "Figure_S2_1.png"), plot=gp, dpi=350, width=2.5, height=4.4)

# Count part
gp <- fit3_pred %>% subset(category=="poi") %>% ggplot(aes(x=vcb))+
  geom_line(aes(y=MED, group=category), color="#D95F02", lwd=2)+
  geom_ribbon(aes(ymax=upr, ymin=lwr), fill="#D95F02", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700))+
  scale_y_continuous(limits=c(0,126))+
  labs(y=expression(paste("Mean number of scale errors (", italic(lambda), ")", sep="")), x="Total vocabulary size")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=0.5),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none") 
print(gp)
ggsave(file = here("Figures", "Figure_S2_2.png"), plot=gp, dpi=350, width=2.5, height=4.4) 

# Overall ZIP
d %>% group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  summarise(Max = max(n), MIN = min(n))   # Max = 3, MIN = 1

gp <- d %>%  group_by(vcb) %>% count(se_sum) %>% ungroup() %>% 
  ggplot(aes(x=vcb))+
  geom_point(aes(y=se_sum, size=n), shape=16, color="gray60", alpha=0.4)+
  geom_line(data=fit3_pred %>% subset(category=="all"), aes(y=MED, group=category), color="#7671B3", lwd=2)+
  geom_ribbon(data=fit3_pred %>% subset(category=="all"), aes(ymax=upr, ymin=lwr), fill="#7671B3", alpha=0.25)+
  theme_classic()+
  scale_x_continuous(limits=c(0,700))+
  scale_y_continuous(limits=c(0,29))+
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
ggsave(file = here("Figures", "Figure_S2_3.png"), plot=gp, dpi=350, width=3.0, height=4.4)  

# Histogram
fit3_dist <- ggs(fit3.bs) %>% 
  filter(str_detect(Parameter, pattern="_dist")) 

d %>% ggplot()+
  geom_histogram(aes(x=se_sum, y=..density..), bins=12, color="white", fill="red", alpha=0.3)+
  geom_histogram(data=fit3_dist, aes(x=value, y=..density..), bins=12, color="white", fill="blue", alpha=0.3)+
  theme_bw()

