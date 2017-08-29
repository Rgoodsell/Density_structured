rm(list=ls())
library(dplyr)
library(rstan)



###### Simulated data function #######
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
}


 
######### Params for simulated data ###############

nsim <- 2e5 # Number of observations
ncat <- 5 #number of simulated categories
pinit <- c(1,0,0,0,0) #initial state probabilities for simulated data - all start as absent here...
alpha <- t(matrix(seq(0,4,by=1),nrow=ncat,ncol=ncat,byrow = T)) # Alpha with 0 constraint
beta <- array(seq(1,10,by=1),c(5,5,2)) # beta
beta[1,,] <- 0 #identifiability constraint
nfield <- 100 # Number of fields
vfield.source <- 1:5 # Field level variance dependent on source state


###### coded explicitly - P does not sum to 1 ###### 
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

#### Using Categorical_logit sampling statement #####
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
# 
# 
# 
# alpha = append_row(alpha_zero,alpha_raw);
# gamma = append_row(gamma_zeroes,gamma_raw);
# 
# 
# 
# 
# 
# 
# }
# 
# model {
# vfield ~ normal(0, 10);
# 
# for(k in 1:K-1){
# alpha_raw[k] ~ normal(0,10);
# for(f in 1:Fi){
# gamma_raw[k,f] ~ normal(0,vfield);
# }}
# 
# 
# 
# 
# 
# for(n in 1:N){
# y[n] ~ categorical_logit(alpha + gamma[,field[n]]);//categorical distribution
# }
# } "

mydata <- simulate.random(n = nsim, ncat = ncat, pstart = pinit, alpha = alpha, nfield = nfield, vfield.source = vfield.source)

# set up for model
y = mydata$ytplus # Source state (all our data = 1 in this case)
N = length(mydata$ytplus) # Number of observations
Fi = nfield # Number of fields 
field = mydata[,3]



field.onlyMod <- stan(model_code =  alpha.gamma.stan, 
                      data=list(N=N,K=5, y = mydata$ytplus, Fi=Fi, field=field), warmup = 1000, iter=2000, 
                      chains = 2, cores=2, control = list(adapt_delta = 0.99,max_treedepth=15))      



traceplot(field.onlyMod,pars=c("alpha","vfield"))
plot(field.onlyMod,pars=c("alpha","gamma"))
pairs(field.onlyMod,pars=c("alpha","vfield"))
summary(field.onlyMod,pars=c("alpha","vfield","gamma"))


saveRDS(field.onlyMod,file = "data/bop14rg/Iceberg_scripts/results/stan_baseline_raneff")
