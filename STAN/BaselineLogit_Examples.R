####### Stan Multinoms ######
## Script contains two stan models
# 1) First model runs a multinomial logistic regression with an intercept and two covariates 
# 2) Second model has an additional random intercept at the field level 

rm(list=ls())
library(dplyr)
library(rstan)

setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/Code/Models/STAN/")



###### Simulated data functions #######

simulate.densitystructure<-function(n,ncat,pstart,alpha,beta){
  #simulate a data set with state yt sampled from a multinomial, and state yt+1 sampled from a baseline-category logit model conditional on yt, with two explanatory variables (whose values are generated from standard normal distributions)
  #n is number of pairs of observations
  #ncat is number of categories
  #pstart is a vector of ncat elements: the parameters of the multinomial distribution for state at time t
  #alpha is an ncat x ncat array, where element alpha_ij is the intercept for category i next year, conditional on category j this year. Identifiability constraint: alpha_0j should equal zero for all j
  #beta is an ncat x ncat x 2 array, where element beta_ijk is the effect of explanatory variable k on the linear predictor for category i next year, conditional on category j this year. Identifiability constraint: beta_0jk should equal 0 for all j and k
  #
  #returns a data frame containing yt,ytplus (simulated density categories at times t,t+1),scalefetch,scalesst (the simulated explanatory variables, sampled from N(0,1),fetch,sst (identical to scalefetch and scalesst: in here only because in the real data, the unscaled explanatory variables don't have standard deviation 1))
  #yt is state this year (sampled from multinomial with probabilities pstart)
  #ytplus is state next year (generated from baseline-category logit model conditional on yt)
  
  yt <- sample(c(1:ncat),replace=T,prob=pstart,size=n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x <- cbind(x1,x2)
  ytplus <- rep(0,n)
  for(j in 1:ncat){#source state
    phi <- array(0,c(n,ncat))
    for(i in 1:ncat){#destination state
      #construct linear predictor for the jth source state and ith destination state
      phi[,i] <- exp(alpha[i,j]+x%*%beta[i,j,])
    }
    p <- phi/rowSums(phi)
    for(l in 1:n){
      if(yt[l]==j){
        ytplus[l] <- sample(c(1:ncat),prob=p[l,],size=1)#sample from categorical distribution with probabilities from baseline-category logit model conditional on current state
      }
    }
  }
  return(data.frame(yt=yt,ytplus=ytplus,scalefetch=x1,scalesst=x2,fetch=x1,sst=x2))
} # Simulates data with non-uniform intercept and continuous  plot level covariate data
simulate.random<- function(n,ncat,nfield,pstart,alpha,vfield.source){
  #simulate a data set with state yt sampled from a multinomial, 
  #and state yt+1 sampled from a baseline-category logit model conditional on yt, 
  #n is number of pairs of observations
  
  #ncat is number of categories
  #pstart is a vector of ncat elements: the parameters of the multinomial distribution for state at time t
  
  #alpha is an ncat x ncat array, where element alpha_ij is the intercept for category i next year,
  #conditional on category j this year.
  #Identifiability constraint: alpha_0j should equal zero for all j
  
  
  #returns a data frame containing yt,ytplus (simulated density categories at times t,t+1),
  #yt is state this year (sampled from multinomial with probabilities pstart)
  #ytplus is state next year (generated from baseline-category logit model conditional on yt)
  
  # builds matrix for yt and field identity and empty matrix for ytplus
  yt <- matrix(rep(0,n), ncol=2, nrow=n)
  yt[,1] <- sample(c(1:ncat),replace=T,prob=pstart,size=n)
  yt[,2] <- sample(c(1:nfield), replace=T,size=n)
  ytplus <- matrix(rep(0,n),ncol=2,nrow=n)
  
  
  gamma <- array(0,dim = c(ncat,ncat,nfield))
  phi <- array(0,dim = c(ncat,ncat,nfield))
  p <- array(0,dim=c(ncat,ncat,nfield))
  
  # generates field effects (gamma), from vfield dependent on each source state
  
  
  for (li in 1:ncat){
    gamma[,li,]<- rnorm(ncat*nfield, 0, sqrt(vfield.source[li]))
  }
  
  gamma[1,,] <- 0
  
  
  
  print(gamma)
  
  print("vfield estimates after setting gamma[1,,]==0  constraint?")
  for (i in 1:ncat){print(var(c(gamma[-1,i,])))}
  
  
  
  for(f in 1:nfield){#field
    for(j in 1:ncat){#source state
      
      for(i in 1:ncat){#destination state
        #construct linear predictor for the jth source state and ith destination state
        phi[i,j,f] <- exp(alpha[i,j]+gamma[i,j,f])
      }}
    
    p[,,f] <- t(t(phi[,,f])/colSums(phi[,,f]))
  }
  
  
  
  
  #sample from categorical distribution with probabilities from baseline-category logit model conditional on current state and field
  for (f in 1:nfield){
    for (j in 1:ncat){# source cat
      for(l in 1:n){
        if(yt[l,1]==j & yt[l,2]==f){
          ytplus[l,1]<- sample(c(1:ncat),prob=p[,j,f],size=1)
          ytplus[l,2] <- f # add in field identity
        }
        
      }
    }}
  return(data.frame(yt=yt[,1],ytplus=ytplus[,1],field=yt[,2]))
}# Simulates data with intercept and random field level intercepts
simulate.random.coeff<- function(n,ncat,nfield,pstart,alpha,beta,vfield.source){
  #simulate a data set with state yt sampled from a multinomial, 
  #and state yt+1 sampled from a baseline-category logit model conditional on yt, 
  #n is number of pairs of observations
  
  #ncat is number of categories
  #pstart is a vector of ncat elements: the parameters of the multinomial distribution for state at time t
  
  #alpha is an ncat x ncat array, where element alpha_ij is the intercept for category i next year,
  #conditional on category j this year.
  #Identifiability constraint: alpha_0j should equal zero for all j
  
  
  #returns a data frame containing yt,ytplus (simulated density categories at times t,t+1),
  #yt is state this year (sampled from multinomial with probabilities pstart)
  #ytplus is state next year (generated from baseline-category logit model conditional on yt)
  
  yt <- matrix(rep(0,n), ncol=4, nrow=n)
  yt[,1] <- sample(c(1:ncat),replace=T,prob=pstart,size=n)
  yt[,2] <- sample(c(1:nfield), replace=T,size=n)
  
  
  # assigns rotational covariates depening on whether field ID is odd or even ¯\_(ツ)_/¯
  for (i in 1:n){
    if(yt[i,2] %%  2 == 0)
    {yt[i,3] <- 1}
    else(yt[i,4] <- 1)}
  
  
  
  ytplus <- matrix(rep(0,n),ncol=4,nrow=n)
  
  # matrix of explanatory variables by field
  x1 <-  rnorm(n)
  x2 <-  rnorm(n)
  x <- cbind(x1,x2)
  
  ytplus <- matrix(rep(0,n),ncol=4,nrow=n)
  
  gamma <- array(0,dim = c(ncat,ncat,nfield))
  phi <- array(0,dim = c(ncat,ncat,nfield))
  p <- array(0,dim=c(ncat,ncat,nfield))
  
  # generates field effects (gamma), from vfield dependent on each source state
  
  for (li in 1:ncat){
    gamma[,li,]<- rnorm(ncat*nfield, 0, sqrt(vfield.source[li]))
  }
  
  # set constraint
  gamma[1,,] <- 0
  
  #print("vfield estimates after setting gamma[1,,]==0  constraint?")
  #for (i in 1:ncat){print(var(c(gamma[,i,])))}
  
  
  
  for(f in 1:nfield){ # field
    for(j in 1:ncat){ # source state 
      for(i in 1:ncat){ # destination state
        #construct linear predictor for jth source state, ith destination state and fth field, for the corresponding coefficients in field f
        phi[i,j,f] <- exp(alpha[i,j] + x%*%beta[i,j,] + gamma[i,j,f])
        
      }
    }
    # calculate transition probabilities 
    p[,,f] <- t(t(phi[,,f])/colSums(phi[,,f]))
    
  }
  
  
  
  
  
  
  
  
  #sample from categorical distribution with probabilities from baseline-category logit model conditional on current state and field
  for (f in 1:nfield){
    for (j in 1:ncat){# source cat
      for(l in 1:n){
        if(yt[l,1]==j & yt[l,2]==f){
          ytplus[l,1]<- sample(c(1:ncat),prob=p[,j,f],size=1)
          ytplus[l,2] <- f # add in field identity
          ytplus[l,3] <- yt[l,3]
          ytplus[l,4] <- yt[l,4]
        }
        
      }
    }}
  return(data.frame(yt=yt[,1],ytplus=ytplus[,1],field=yt[,2],rt1=ytplus[,3],rt2=ytplus[,4]))
} # Simulate data with categorical field level covariates and random intercepts at the field level 

######### Params for simulated data ###############

nsim <- 1e3 # Number of observations
ncat <- 5 #number of simulated categories
pinit <- c(1,0,0,0,0) #initial state probabilities for simulated data - all start as absent here...
alpha <- t(matrix(seq(0,4,by=1),nrow=ncat,ncol=ncat,byrow = T)) # Alpha with 0 constraint
beta <- array(seq(1,10,by=1),c(5,5,2)) # beta
beta[1,,] <- 0 #identifiability constraint
nfield <- 10 # Number of fields
vfield.source <- 1:5 # Field level variance dependent on source state - check these after simulating data, likely to change after setting constraints on random effect


####### Intercept + gradient stan model ######
##############################################

# Model code
alpha.beta.stan <-  "data {
int<lower=2> K; // No. of categories
int<lower=0> N; // No. of  observations
int<lower=1> D; // No. of explanatory vars
int<lower=1,upper=K> y[N]; // State outcome
vector[D] x[N]; //  Explanatory variable matrix
}

transformed data {
row_vector[D] beta_zeroes;
int<lower=0> alpha_zero;
beta_zeroes = rep_row_vector(0,D); 
alpha_zero = 0;

}
parameters {
vector[K-1] alpha_raw;
matrix[K-1,D] beta_raw;

}

transformed parameters{
matrix[K,D] beta;
vector[K] alpha;

alpha = append_row(alpha_zero,alpha_raw);
beta = append_row(beta_zeroes, beta_raw);

}

model {

for(k in 1:K-1){
alpha_raw[k] ~ normal(0,10);
beta_raw[k] ~ normal(0,10);
}

for(n in 1:N){

y[n] ~ categorical_logit(alpha + (beta * x[n])); 

}

} "

# simulate data and check
mydata<- simulate.densitystructure(n = nsim ,ncat = ncat,pstart = pinit,alpha = alpha,beta=beta)
str(mydata)

# set up data for model
y = mydata$ytplus # Source state (all our data = 1 in this case)
x = as.matrix(mydata[3:4]) # Explanatory variables
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 



# Run the model
igMod <- stan(model_code =  alpha.beta.stan, 
                     data=list(N=N,K=5,D=D,x=x, y = y), warmup = 500, iter=2000, 
                     chains = 2, cores=2, control = list(adapt_delta = 0.99))


# Check the outputs
traceplot(igMod,pars=c("alpha","beta"))
plot(igMod,pars=c("alpha","beta"))
pairs(igMod,pars=c("alpha","beta"))
summary(igMod)



########## Intercept and field level intercept model ##########
################################################


# simulate data 
mydata <- simulate.random(n = nsim, ncat = ncat, pstart = pinit, alpha = alpha, nfield = nfield, vfield.source = vfield.source)

# set up for model
y = mydata$ytplus # Source state (all our data = 1 in this case)
N = length(mydata$ytplus) # Number of observations
Fi = nfield # Number of fields 
field = mydata[,3]



###### baseline logit model coded explicitly ######

alpha.gamma.stan <-  "data {
int<lower=2> K; // No. of categories
int<lower=0> N; // No. of  observations
int<lower=1> Fi; // No. of fields (random intercepts)
int<lower=1,upper=K> y[N]; // State outcome
int<lower=1> field[N]; // map fields to observations
}

transformed data {
int<lower=0> alpha_zero;
row_vector[Fi] gamma_zeroes;
row_vector[N] phi_1;

gamma_zeroes = rep_row_vector(0,Fi);
alpha_zero = 0;
phi_1 = rep_row_vector(1,N);

}
parameters {
vector[K-1] alpha_raw;
matrix[K-1,Fi] gamma_raw;
real<lower=0> vfield;
matrix[K-1,N] phi_raw;

}

transformed parameters{
vector[K] alpha;
matrix[K,Fi] gamma;
matrix[K,N] phi;
matrix[K,N] p;

alpha = append_row(alpha_zero,alpha_raw);
phi = append_row(phi_1, phi_raw);
gamma = append_row(gamma_zeroes,gamma_raw);

for(n in 1:N){
phi[,n] = exp(alpha + gamma[,field[n]]);
p[,n] = phi[,n]/sum(phi[,n]);//baseline-category logit

}

}

model {
vfield ~ normal(0, 10);

for(k in 1:K-1){
alpha_raw[k] ~ normal(0,10);
for(f in 1:Fi){
gamma_raw[k,f] ~ normal(0,vfield);
}}

for(n in 1:N){
y[n] ~ categorical(p[,n]);//categorical distribution
}
} "



### Using categorical_logit sampling statement ### 
# alpha.gamma.stan <-   "data {
# int<lower=2> K; // No. of categories
# int<lower=0> N; // No. of  observations
# int<lower=1> Fi; // No. of fields (random intercepts)
# int<lower=1,upper=K> y[N]; // State outcome
# int<lower=1> field[N]; // map fields to observations
# }
# 
# transformed data {
# int<lower=0> alpha_zero;
# row_vector[Fi] gamma_zeroes;
# 
# alpha_zero = 0;
# gamma_zeroes = rep_row_vector(0,Fi);
# 
# }
# parameters {
# vector[K-1] alpha_raw;
# matrix[K-1,Fi] gamma_raw;
# real<lower=0> vfield;
# 
# }
# 
# transformed parameters{
# vector[K] alpha;
# matrix[K,Fi] gamma;
# alpha = append_row(alpha_zero,alpha_raw);
# gamma = append_row(gamma_zeroes,gamma_raw);
# }
# 
# model {
# vfield ~ normal(0, 10);

# for(k in 1:K-1){
# alpha_raw[k] ~ normal(0,10);
# for(f in 1:Fi){
# gamma_raw[k,f] ~ normal(0,vfield);
# }}
# 
# for(n in 1:N){
# y[n] ~ categorical_logit(alpha + gamma[,field[n]]);//categorical distribution
# }
# } "


field.onlyMod <- stan(model_code =  alpha.gamma.stan, 
                      data=list(N=N,K=5, y = mydata$ytplus, Fi=Fi, field=field), warmup = 1000, iter=2000, 
                      chains = 2, cores=2, control = list(adapt_delta = 0.99))      



traceplot(field.onlyMod,pars=c("alpha","vfield"))
plot(field.onlyMod,pars=c("alpha","vfield"))
pairs(field.onlyMod,pars=c("alpha","vfield"))
summary(field.onlyMod,pars=c("alpha","vfield","gamma"))



#### Random effect and explanatory vars stan model  #####
#############################################

# Model code
alpha.beta.gamma.stan <-  "data {
int<lower=2> K; // No. of categories
int<lower=0> N; // No. of  observations
int<lower=1> D; // No. of explanatory vars
int<lower=1> Fi; // No. of fields (random intercepts)
int<lower=1,upper=K> y[N]; // State outcome
vector[D] x[N]; //  Explanatory variable matrix
int<lower=1> field[N]; // map fields to observations
}

transformed data {
row_vector[D] beta_zeroes;
int<lower=0> alpha_zero;
row_vector[Fi] gamma_zeroes;
beta_zeroes = rep_row_vector(0,D); 
alpha_zero = 0;
gamma_zeroes = rep_row_vector(0,Fi);

}
parameters {
vector[K-1] alpha_raw;
matrix[K-1,D] beta_raw;
matrix[K-1,Fi] gamma_raw;
real<lower=0> vfield;

}

transformed parameters{
matrix[K,D] beta;
vector[K] alpha;
matrix[K,Fi] gamma;



alpha = append_row(alpha_zero,alpha_raw);
beta = append_row(beta_zeroes, beta_raw);
gamma = append_row(gamma_zeroes,gamma_raw);

}

model {
vfield ~ normal(0, 5);

for(k in 1:K-1){
alpha_raw[k] ~ normal(0,10);
beta_raw[k,] ~ normal(0,10);
gamma_raw[k,] ~ normal(0,vfield);
}



for(n in 1:N){
y[n] ~ categorical_logit((alpha + gamma[,field[n]]) + (beta * x[n])) ; 

}

} "

# simulate data and check
mydata<- simulate.random.coeff(n = nsim ,ncat = ncat,pstart = pinit,alpha = alpha,beta=beta,nfield = nfield,vfield.source = vfield.source)

# set up for model
y = mydata$ytplus # Source state (all our data = 1 in this case)
x = as.matrix(mydata[3:4]) # Explanatory variables
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 
Fi = nfield # Number of fields 
field = mydata[,3]

# Run model 
field.interceptMod <- stan(model_code =  alpha.beta.gamma.stan , 
                           data=list(N=N,K=5,D=D,x=x, y = y, Fi=Fi, field=field), warmup = 500, iter=2000, 
                           chains = 2, cores=2, control = list(adapt_delta = 0.99))


# Check outputs 
traceplot(field.interceptMod,pars=c("alpha","beta","vfield"))
plot(field.interceptMod,pars=c("alpha","beta","vfield"))
pairs(field.interceptMod,pars=c("alpha","beta","vfield"))
summary(field.interceptMod)

                           
                           

