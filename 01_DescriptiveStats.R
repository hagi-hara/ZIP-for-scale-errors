rm(list=ls())

R.version           # version 4.2.2
library(tidyverse)  # version 2.0.0

# read raw data
d <- read.csv("sedata.csv", header=TRUE)

# Histogram of SE
d %>% subset(setting=="inlab") %>% 
  ggplot(aes(x=se_sum, y=..density..))+
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
d.par <- d %>% select(gender, parid, groupid, setting, country) %>% unique()
d.par %>% nrow()
d.par$gender %>% table()
d.par$country %>% table()

d %>% group_by(setting) %>% summarize(N = n())
d %>% group_by(setting, country) %>% summarize(N = n())
d.par %>% group_by(setting, country) %>% summarize(N = n())

d %>% summarize(N = n(),
                MeanAge  = mean(age),
                SDAge    = sd(age),
                MinAge   = min(age),
                MaxAge   = max(age),
                nGirl    = sum(gender =="f", na.rm=TRUE),
                propGirl = nGirl/N,
                nSE      = sum(se_occ=="1"),
                PropSE   = nSE/N,
                MeanSE   = mean(se_sum),
                SDSE     = sd(se_sum),
                MinSE    = min(se_sum),
                MaxSE    = max(se_sum))
  
d %>% group_by(groupid, setting) %>% 
  summarize(N = n(),
            MeanAge  = mean(age, na.rm=TRUE),
            SDAge    = sd(age, na.rm=TRUE),
            MinAge   = min(age, na.rm=TRUE),
            MaxAge   = max(age, na.rm=TRUE),
            nGirl    = sum(gender=="f", na.rm=TRUE),
            propGirl = nGirl/N,
            nSE      = sum(se_occ=="1"),
            PropSE   = nSE/N,
            MeanSE   = mean(se_sum),
            SDSE     = sd(se_sum),
            MinSE    = min(se_sum),
            MaxSE    = max(se_sum))

d %>% subset(setting=="inlab") %>% group_by(gender, session) %>% 
  summarize(N = n(),
            MeanAge  = mean(age, na.rm=TRUE),
            SDAge    = sd(age, na.rm=TRUE),
            MinAge   = min(age, na.rm=TRUE),
            MaxAge   = max(age, na.rm=TRUE),
            nSE      = sum(se_occ=="1"),
            PropSE   = nSE/N,
            MeanSE   = mean(se_sum),
            SDSE     = sd(se_sum),
            MinSE    = min(se_sum),
            MaxSE    = max(se_sum))

# Data with classroom settings
d.cl <- d %>% subset(setting=="classroom")

d.cl %>% select(parid, groupid, gender) %>% unique() %>% 
  group_by(groupid) %>% 
  summarize(N = n(),
            nGirl    = sum(gender=="f", na.rm=TRUE),
            propGirl = nGirl/N)

d.cl %>% select(parid, groupid, gender) %>% unique() %>% 
  subset(groupid!="R2009") %>% 
  summarize(N = n(),
            nGirl    = sum(gender=="f", na.rm=TRUE),
            propGirl = nGirl/N)

d.cl %>% group_by(parid, groupid) %>% summarize(NSes = n()) %>% 
  group_by(groupid) %>% 
  summarize(N       = n(),
            MeanSes = mean(NSes, na.rm=TRUE),
            SDSes   = sd(NSes, na.rm=TRUE), 
            MinSes  = min(NSes, na.rm=TRUE),
            MaxSes  = max(NSes, na.rm=TRUE))


# read raw data including vocabulary size 
d2 <- d %>% subset(vcball!="NA") 

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
                 nGirl    = sum(gender =="f", na.rm=TRUE),
                 propGirl = nGirl/N,
                 nSE      = sum(se_occ=="1"),
                 PropSE   = nSE/N,
                 MeanSE   = mean(se_sum),
                 SDSE     = sd(se_sum),
                 MinSE    = min(se_sum),
                 MaxSE    = max(se_sum),
                 MeanVcv  = mean(vcball),
                 SDVcb    = sd(vcball),
                 MinVcb   = min(vcball),
                 MaxVcb   = max(vcball))

d2 %>% group_by(groupid) %>% 
  summarize(N = n(),
            MeanAge  = mean(age),
            SDAge    = sd(age),
            MinAge   = min(age),
            MaxAge   = max(age),
            nGirl    = sum(gender =="f", na.rm=TRUE),
            propGirl = nGirl/N,
            MeanVcv  = mean(vcball),
            SDVcb    = sd(vcball),
            MinVcb   = min(vcball),
            MaxVcb   = max(vcball))

d2 %>% subset(noun!="NA") %>% 
  summarize(N = n(),
            MeanNoun = mean(noun, na.rm=TRUE),
            SDNoun   = sd(noun, na.rm=TRUE),
            MinNoun  = min(noun, na.rm=TRUE),
            MaxNoun  = max(noun, na.rm=TRUE),
            MeanVerb = mean(verb, na.rm=TRUE),
            SDVerb   = sd(verb, na.rm=TRUE),
            MinVerb  = min(verb, na.rm=TRUE),
            MaxVerb  = max(verb, na.rm=TRUE),
            MeanAdj  = mean(adj, na.rm=TRUE),
            SDAdj    = sd(adj, na.rm=TRUE),
            MinAdj   = min(adj, na.rm=TRUE),
            MaxAdj   = max(adj, na.rm=TRUE))

d2 %>% subset(noun!="NA") %>% group_by(groupid) %>% 
  summarize(N = n(),
            MeanNoun = mean(noun, na.rm=TRUE),
            SDNoun   = sd(noun, na.rm=TRUE),
            MinNoun  = min(noun, na.rm=TRUE),
            MaxNoun  = max(noun, na.rm=TRUE),
            MeanVerb = mean(verb, na.rm=TRUE),
            SDVerb   = sd(verb, na.rm=TRUE),
            MinVerb  = min(verb, na.rm=TRUE),
            MaxVerb  = max(verb, na.rm=TRUE),
            MeanAdj  = mean(adj, na.rm=TRUE),
            SDAdj    = sd(adj, na.rm=TRUE),
            MinAdj   = min(adj, na.rm=TRUE),
            MaxAdj   = max(adj, na.rm=TRUE))

d2 %>% subset(noun_sp!="NA") %>% 
  summarize(N = n(),
            MeanNoun = mean(noun_sp, na.rm=TRUE),
            SDNoun   = sd(noun_sp, na.rm=TRUE),
            MinNoun  = min(noun_sp, na.rm=TRUE),
            MaxNoun  = max(noun_sp, na.rm=TRUE),
            MeanVerb = mean(verb_sp, na.rm=TRUE),
            SDVerb   = sd(verb_sp, na.rm=TRUE),
            MinVerb  = min(verb_sp, na.rm=TRUE),
            MaxVerb  = max(verb_sp, na.rm=TRUE),
            MeanAdj  = mean(adj_sp, na.rm=TRUE),
            SDAdj    = sd(adj_sp, na.rm=TRUE),
            MinAdj   = min(adj_sp, na.rm=TRUE),
            MaxAdj   = max(adj_sp, na.rm=TRUE))

d2 %>% subset(noun_sp!="NA") %>% group_by(groupid) %>% 
  summarize(N = n(),
            MeanNoun = mean(noun_sp, na.rm=TRUE),
            SDNoun   = sd(noun_sp, na.rm=TRUE),
            MinNoun  = min(noun_sp, na.rm=TRUE),
            MaxNoun  = max(noun_sp, na.rm=TRUE),
            MeanVerb = mean(verb_sp, na.rm=TRUE),
            SDVerb   = sd(verb_sp, na.rm=TRUE),
            MinVerb  = min(verb_sp, na.rm=TRUE),
            MaxVerb  = max(verb_sp, na.rm=TRUE),
            MeanAdj  = mean(adj_sp, na.rm=TRUE),
            SDAdj    = sd(adj_sp, na.rm=TRUE),
            MinAdj   = min(adj_sp, na.rm=TRUE),
            MaxAdj   = max(adj_sp, na.rm=TRUE))


d2 %>% subset(noun!="NA") %>% 
  ggplot(aes(x=verb, y=se_occ))+
  geom_point()+
  stat_smooth(method=lm, formula=y~x+I(x^2), se=TRUE, lwd=2)
