##### Ordered category Logit in Stan #####

rm(list=ls())

setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/Code/Models/STAN/")


library(rstan)
library(dplyr)


###################

##  Data Setup  ##

###################

# Read in data, dop nas & unused levels
rawData <- read.csv("rotation_2007_2008_2009_2010.csv") %>% na.omit(droplevels(.))

str(rawData)


# Select fields with lots of high DS
highStates<- droplevels(mydata[rawData$Year1 > 3,])
highFields <- levels(highStates$field) # Select random intercept level here, fields or farms
highData <- droplevels(rawData[rawData$field %in% highFields,])


# setup for jags models 
mydata <-  droplevels(highData) %>%
  mutate(yt = factor(Year1 +1), 
         ytplus = Year2 +1) %>% 
  select (yt, ytplus,field) 
  
str(mydata)

########################

## Fixed Effects Model ##

#########################

ordered.logisitc.stan <-  "data {
  int<lower=2> K; // Categories
  int<lower=0> N; // No. Observations
  int<lower=1> D; // Explanatory Variables (source state)
  int<lower=1,upper=K> y[N]; // Outcome
  row_vector[D] x[N]; // Explanatory variable design matrix
}

transformed data {
int<lower=0> beta_zero;
beta_zero = 0;
}

parameters {
  vector[D-1] beta_raw; // Source state effect
  ordered[K-1] c; // Cutpoints


}

transformed parameters{ // K-1 parameterisation
 vector[D] beta;
  beta = append_row(beta_zero,beta_raw);
}
model {
  for (n in 1:N)
    y[n] ~ ordered_logistic(x[n] * beta, c);
}"



y = mydata$ytplus # Outcome
x <- data.frame(model.matrix( ~yt -1, data = mydata)) # Design matrix for source state
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 

# Run the model
ordered.logitMod <- stan(model_code =  ordered.logisitc.stan, 
              data=list(N=N,K=5,D=D,x=x, y = y), warmup = 500, iter=1000, 
              chains = 2, cores=2, control = list(adapt_delta = 0.99))


# Check the outputs
traceplot(ordered.logitMod,pars=c("beta","c"))
plot(ordered.logitMod,pars=c("beta","c"))
pairs(ordered.logitMod,pars=c("beta","c"))
summary(ordered.logitMod)


################################

## Random field effects Model ##

################################

ordered.raneff.multi.stan <-  "data {
int<lower=2> K;
int<lower=0> N;
int<lower=1> D;
int<lower=1> Fi; // Number of groups (fields or farms)
int<lower=1,upper=K> y[N];
row_vector[D] x[N];
int<lower=1,upper=Fi> field[N]; // Indexing variable
}

transformed data {
int<lower=0> beta_zero;
row_vector[Fi] gamma_zeroes;

gamma_zeroes = rep_row_vector(0,Fi);
beta_zero = 0;
}

parameters {
vector[D-1] beta_raw;
matrix[D-1,Fi] gamma_raw; // Random intercept
vector<lower=0>[D-1] sigmaGamma;  // Varianace for random intercept
ordered[K-1] c;


}

transformed parameters{
vector[D] beta;
matrix[D,Fi] gamma; // Random field effect
matrix[D,Fi] eta;
beta = append_row(beta_zero,beta_raw);
gamma = append_row(gamma_zeroes,gamma_raw);

for(d in 1:D){
for(f in 1:Fi){
eta[d,f] = beta[d] + gamma[d, f]; // constructing eta parameter here lets things run smoother for some reason....
}}
}
model {




for(d in 1:D-1){
sigmaGamma[d] ~ cauchy(0,5);
beta_raw[d] ~ normal(0,5);
for(f in 1:Fi){
gamma_raw[d,f] ~ normal(0,sigmaGamma);
}}


for (n in 1:N){
y[n] ~ ordered_logistic(x[n] * eta[,field[n]], c);
}
}"





y = mydata$ytplus # Source state (all our data = 1 in this case)
x <- data.frame(model.matrix( ~yt -1, data = mydata))
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 
Fi = length(unique(mydata$field))
field = as.numeric(droplevels(mydata$field))

ordered.raneffMod <- stan(model_code =  ordered.raneff.stan, 
                         data=list(N=N,K=5,D=D,x=x, y = y,Fi=Fi,field=field), warmup = 1000, iter=2000, 
                         chains = 2, cores=2, control = list(adapt_delta = 0.99,max_treedepth = 15))


traceplot(ordered.raneffMod,pars=c("beta"))
pairs(ordered.raneffMod,pars=c("beta"))


# Check a couple of field level intercepts
traceplot(ordered.raneffMod,pars=c("gamma[2,5]"))
plot(ordered.raneffMod,pars=c("gamma[2,20"))





