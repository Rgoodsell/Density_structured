##### Ordered category Logit in Stan #####

rm(list=ls())

setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/")


library(rstan)
library(dplyr)



###################
###################
# read in
mydata <- read.csv("data/rotation_2007_2008_2009_2010.csv")

str(mydata)

# select 4 most common rotations
mydata <- na.omit (droplevels(mydata)) 

str(mydata)


# setup for jags models 
mydata <- droplevels(mydata) %>% 
  
  mutate(yt = Year1 +1 , 
         ytplus = Year2 +1) %>% 
  
  select (yt, ytplus,field) %>% 
  sample_n(10000)

mydata$yt <- factor(mydata$yt)

######################
#######################



ordered.logisitc.stan <-  "data {
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  row_vector[D] x[N];
}

transformed data {
int<lower=0> beta_zero;
beta_zero = 100;
}

parameters {
  vector[D-1] beta_raw;
  ordered[K-1] c;


}

transformed parameters{
 vector[D] beta;
  beta = append_row(beta_zero,beta_raw);
}
model {
for(d in 1:D-1)
beta_raw[d] ~ normal(0,100);

  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, c);
}"



y = mydata$ytplus # Source state (all our data = 1 in this case)
x <- data.frame(model.matrix( ~yt -1, data = mydata))
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 

# Run the model
ordered.logitMod <- stan(model_code =  ordered.logisitc.stan, 
              data=list(N=N,K=5,D=D,x=x, y = y), warmup = 1000, iter=1000, 
              chains = 2, cores=2, control = list(adapt_delta = 0.99))


# Check the outputs
traceplot(ordered.logitMod,pars=c("beta"))
plot(ordered.logitMod,pars=c("beta"))
pairs(ordered.logitMod,pars=c("beta"))
summary(ordered.logitMod)
# 
# 
# ordered.logitMod.e <-  extract(ordered.logitMod)
# 
# #cutpoints_B1_0<- apply(ordered.logitMod.e$c,2,mean)
# #Beta_B1_0 <- apply(ordered.logitMod.e$beta,2,mean)
# 
# #cutpoints_B1_1 <- apply(ordered.logitMod.e$c,2,mean)
# #Beta_B1_1 <- apply(ordered.logitMod.e$beta,2,mean)
# 
# #cutpoints_B1_10 <- apply(ordered.logitMod.e$c,2,mean)
# #Beta_B1_10 <- apply(ordered.logitMod.e$beta,2,mean)
# 
# cutpoints_B1_100 <- apply(ordered.logitMod.e$c,2,mean)
# Beta_B1_100 <- apply(ordered.logitMod.e$beta,2,mean)
# 
# 
# 
# betas <- c(Beta_B1_0[-1],Beta_B1_1[-1],Beta_B1_10[-1],Beta_B1_100[-1])
# cutpoints <- c(cutpoints_B1_0,cutpoints_B1_1+1,cutpoints_B1_10,cutpoints_B1_100)
# 
# 
# 
# betas <- data.frame(betas,B1=c(0,0,0,0,1,1,1,1,10,10,10,10,100,100,100,100),cp=rep(c(2,3,4,5),4))
# cutpoints <- data.frame(cutpoints, B1=c(0,0,0,0,1,1,1,1,10,10,10,10,100,100,100,100),cp=rep(c(1,2,3,4),4))
# 
# 
# 
# 
# 
# 
# ggplot(cutpoints,aes(B1,cutpoints))+
#   geom_point(aes(colour=factor(cp),size=1))
# 
# ggplot(betas,aes(B1,betas))+
#   geom_point(aes(colour=factor(cp),size=1))
# 





ordered.raneff.stan <-  "data {
int<lower=2> K;
int<lower=0> N;
int<lower=1> D;
int<lower=1> Fi;
int<lower=1,upper=K> y[N];
row_vector[D] x[N];
int<lower=1> field[N];
}

transformed data {
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
vector[Fi] gamma;
ordered[K-1] c;


}

transformed parameters{
vector[D] beta;
beta = append_row(beta_zero,beta_raw);
}
model {
for(d in 1:D-1)
beta_raw[d] ~ normal(0,5);

for (n in 1:N)
y[n] ~ ordered_logistic(x[n] * beta + gamma[field[n]], c);
}"




y = mydata$ytplus # Source state (all our data = 1 in this case)
x <- data.frame(model.matrix( ~yt -1, data = mydata))
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 
Fi = length(unique(mydata$field))
field = as.numeric(mydata$field)

ordered.ranefftMod <- stan(model_code =  ordered.raneff.stan, 
                         data=list(N=N,K=5,D=D,x=x, y = y,Fi=Fi,field=field), warmup = 1000, iter=1000, 
                         chains = 2, cores=2, control = list(adapt_delta = 0.99))



