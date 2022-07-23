R.version           # version 4.0.4 
library(tidyverse)  # version 1.3.0

# read raw data
d <- read_csv("sedata.csv", col_types="ccccdcdccd")

# Histogram of SE
d %>% ggplot(aes(x=se_sum, y=..density..))+
  geom_histogram(bins=8, color="white")+
  theme_bw()

mean(d$se_sum)
var(d$se_sum)

seed=123
d.sim <- data.frame(se_sum=rpois(100000, lambda=mean(d$se_sum)))
d.sim <- data.frame(se_sum=rpois(100000, lambda=0.8))
d %>% ggplot(aes(x=se_sum, y=..density..))+
  geom_histogram(bins=8, color="white", fill="red", alpha=0.5)+
  geom_histogram(data=d.sim, bins=8, color="white", fill="blue", alpha=0.5)+
  theme_bw()


# Participant characteristics
d %>% summarize(N = n(),
                MeanAge  = mean(age),
                SDAge    = sd(age),
                MinAge   = min(age),
                MaxAge   = max(age),
                nGirl    = sum(gender =="f"),
                propGirl = nGirl/N,
                nSE      = sum(se_occ=="1"),
                PropSE   = nSE/N,
                MeanSE   = mean(se_sum),
                SDSE     = sd(se_sum),
                MinSE    = min(se_sum),
                MaxSE    = max(se_sum))
  
d %>% group_by(groupid) %>% 
  summarize(N = n(),
            MeanAge  = mean(age),
            SDAge    = sd(age),
            MinAge   = min(age),
            MaxAge   = max(age),
            nGirl    = sum(gender =="f"),
            propGirl = nGirl/N,
            nSE      = sum(se_occ=="1"),
            PropSE   = nSE/N,
            MeanSE   = mean(se_sum),
            SDSE     = sd(se_sum),
            MinSE    = min(se_sum),
            MaxSE    = max(se_sum))

# read raw data including vocabulary size 
d2 <- d %>% subset(vcb!="NA") 

nrow(d2)/nrow(d)
  
# Histogram of SE
d2 %>% ggplot(aes(x=se_sum, y=..density..))+
  geom_histogram(bins=8, color="white")+
  theme_bw()

mean(d2$se_sum)
var(d2$se_sum)

seed=123
d.sim <- data.frame(se_sum=rpois(100000, lambda=mean(d$se_sum)))
d.sim <- data.frame(se_sum=rpois(100000, lambda=0.8))
d2 %>% ggplot(aes(x=se_sum, y=..density..))+
  geom_histogram(bins=8, color="white", fill="red", alpha=0.5)+
  geom_histogram(data=d.sim, bins=8, color="white", fill="blue", alpha=0.5)+
  theme_bw()

# Participant characteristics
d2 %>% summarize(N = n(),
                 MeanAge  = mean(age),
                 SDAge    = sd(age),
                 MinAge   = min(age),
                 MaxAge   = max(age),
                 nGirl    = sum(gender =="f"),
                 propGirl = nGirl/N,
                 nSE      = sum(se_occ=="1"),
                 PropSE   = nSE/N,
                 MeanSE   = mean(se_sum),
                 SDSE     = sd(se_sum),
                 MinSE    = min(se_sum),
                 MaxSE    = max(se_sum),
                 MeanVcv  = mean(vcb),
                 SDVcb    = sd(vcb),
                 MinVcb   = min(vcb),
                 MaxVcb   = max(vcb))

d2 %>% group_by(groupid) %>% 
  summarize(N = n(),
            MeanAge  = mean(age),
            SDAge    = sd(age),
            MinAge   = min(age),
            MaxAge   = max(age),
            nGirl    = sum(gender =="f"),
            propGirl = nGirl/N,
            nSE      = sum(se_occ=="1"),
            PropSE   = nSE/N,
            MeanSE   = mean(se_sum),
            SDSE     = sd(se_sum),
            MinSE    = min(se_sum),
            MaxSE    = max(se_sum),
            MeanVcv  = mean(vcb),
            SDVcb    = sd(vcb),
            MinVcb   = min(vcb),
            MaxVcb   = max(vcb))

