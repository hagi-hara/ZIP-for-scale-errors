rm(list=ls())

library(here)       # version 1.0.1
library(splines2)   # version 0.5.0
library(tidyverse)  # version 2.0.0

output_fig_dir <- here("Figures")
if (!dir.exists(output_fig_dir)){
  dir.create(output_fig_dir, recursive = TRUE)
}

# Figure S1

## Simple linear regression
x <- seq(0, 1, 0.01)
set.seed(123)
gp <- data.frame(x = x, y = -6*x+7+rnorm(x,0,1)) %>% 
  ggplot(aes(x=x, y=y))+
  geom_point(size=3, color="gray60")+
  stat_function(fun=function(x) -6*x+7, color="#d11141", size=2.5)+
  scale_y_continuous(limits=c(0,9))+
  theme_classic()+
  labs(y="", x="")+
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(color="black", size=1))
print(gp)
ggsave(file = here("Figures", "Figure_S1_1.png"), plot=gp, dpi=350, width=6.2, height=4.2)    


## Quadratic regression
set.seed(123)
gp <- data.frame(x = x, y = -18*(x-0.5)^2+6+rnorm(x,0,1)) %>% 
  ggplot(aes(x=x, y=y))+
  geom_point(size=3, color="gray60")+
  stat_function(fun=function(x) -18*(x-0.5)^2+6, color="#d11141", size=2.5)+
  scale_y_continuous(limits=c(0,9))+
  theme_classic()+
  labs(y="", x="")+
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(color="black", size=1))
print(gp)
ggsave(file = here("Figures", "Figure_S1_2.png"), plot=gp, dpi=350, width=6.2, height=4.2)    


## B-Spline regression
### base functions
knots <- seq(0, 1, length=3)
knots <- knots[2:(length(knots)-1)]
bsMat <- t(bSpline(x, knots=knots, degree=3, intercept=TRUE))

### Visualization of base functions
gp <- as.data.frame(t(bsMat)) %>% cbind(x) %>% 
  pivot_longer(!`x`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=x, y=value, color=basis_id))+
  geom_line(size=1.2, alpha=0.8)+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray70", size=0.6)+
  scale_y_continuous(limits=c(0,1))+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  labs(y="", x="")+
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(color="black", size=1),
        legend.position="none")
print(gp)
ggsave(file = here("Figures", "Figure_S1_3.png"), plot=gp, dpi=350, width=6.2, height=2.2)

### weights of base functions
w <- c(6, 11, 2, 7, 1)

### Make the figure
set.seed(123)
gp <- as.data.frame(t(w*bsMat)) %>% cbind(x) %>% 
  pivot_longer(!`x`, names_to="basis_id", values_to="value") %>% 
  ggplot(aes(x=x, y=value, color=basis_id))+
  geom_point(data=data.frame(x = x, y = t(w%*%bsMat)+rnorm(x, 0, 0.8)), aes(x=x, y=y), size=3, color="gray60")+
  geom_line(size=1.2, alpha=0.3)+
  geom_vline(data=data.frame(knots), aes(xintercept=knots), linetype="dashed", color="gray70", size=0.6)+
  geom_line(data=data.frame(x = x, y = t(w%*%bsMat)), aes(x=x, y=y), color="#d11141", size=2.5)+
  scale_y_continuous(limits=c(0,10))+
  scale_color_brewer(palette="Dark2")+
  theme_classic()+
  labs(y="", x="")+
  theme(axis.text.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_line(color="black", size=1),
        legend.position="none")
print(gp)
ggsave(file = here("Figures", "Figure_S1_4.png"), plot=gp, dpi=350, width=6.2, height=4.2)



# Figure 1
d.sim <- data.frame(matrix(rep(NA, 3), nrow=1))[numeric(0), ]
colnames(d.sim) <- c("Poisson", "Bernoulli", "ZIP")
set.seed(123)
for(i in 1:10000){
  tmp1 <- rpois(1, 0.9)
  tmp2 <- rbinom(1, 1, 0.7)
  if(tmp2 == 0){tmp3 <- 0} else {tmp3 <- tmp1}
  d.sim[i,] <- c(tmp1, tmp2, tmp3)
}
rm(tmp1, tmp2, tmp3)

d.sim %>% select(ZIP) %>% unique()

d.sim$Poisson %>% mean()    # lambda
d.sim$Bernoulli %>% mean()  # theta

gp <- d.sim %>% mutate(Bernoulli = as.factor(Bernoulli)) %>% 
  ggplot(aes(x=ZIP, fill=Bernoulli))+
  geom_histogram(bins=7, color="white")+
  scale_x_continuous(breaks=0:7)+
  scale_fill_brewer(palette="Dark2")+
  theme_classic()+
  labs(y="Count", x="Outcome")+
  theme(axis.ticks=element_line(color = "black"),
        axis.text.x=element_text(size=12, color="black", angle=0, hjust=1),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14, color="black"),
        strip.text=element_text(size=14, face="bold"),
        legend.position="none")
print(gp)
ggsave(file = here("Figures", "Figure_1.png"), plot=gp, dpi=350, width=4.2, height=4.2)
