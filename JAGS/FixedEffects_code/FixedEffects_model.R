################################################################################################################################################ 
#                                                                                                                                             #
#                                                 ### Density Structured models 4/11/2016 ###                                                 #
#                                                                                                                                             #
#                                                              ## Contents ##                                                                 #
#                                                                                                                                             #
#                                                          # 1) Data & setup                                                                  #
#                                                          # 2) Fixed effects model                                                           #
#                                                                                                                                             #
################################################################################################################################################


rm(list=ls())
setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/")# set path to the distribution folder



################################################################################################################################################ 
#                                                                                                                                             #
#                                                 ### DATA & SETUP ###                                                                        #
#                                                                                                                                             #
################################################################################################################################################



library(nnet) 
library(rjags)
library(mgcv)
library(MASS)
library(dplyr)


# source needed functions
source("simulate.densitystructure.R")
source("chainfuncsmulti.R")
source("construct.P.R")
source("transition.count.R")
source("transition.count.R")
source("Psummary.R")
source("plotP.R")


# I've attached the full dataset and a shorter test dataset that will run much quicker 

# Read in full dataset
#rotation  <-  read.csv(file = "rotation_2007_2008_2009_2010.csv")

# Read in test dataset ( approx. 5000 obvservations selected randomly)
rotation <- na.omit(read.csv(file = "Dylan_test_data.csv"))


# Pretty self explanatory, Year 1 & 2 are density states in respective years, rotation is rotation status (0 = unrotated, 1 = rotated)
# survey is initial survey year(i.e. what year were density states in year1 recorded, year2 states are therefore the subsequent year)
# 
str(rotation)
head(rotation)





################################################################################################################################################ 
#                                                                                                                                             #
#                                                 ### FIXED EFFECTS MODEL ###                                                                 #
#                                                                                                                                             #
################################################################################################################################################



###  Convert col names etc to those that the model recognises ###

 # yt = density state in first year
 # ytplus = density state in second year
 # rotation = rotation status (0 = unrotated, 1 = rotated)

mydata <-  rotation %>% 
  transmute(yt = Year1 +1 , ytplus = Year2 +1, rotation=as.numeric(rotation))
head(mydata)


####################  Test numbers   #################### 

# Matt Spencer mentions that are nowhere near large enough but they give reasonable transition estimates,
# Much larger numbers become unworkable with the number of observations we've got in the full dataset (because of the size of the P array, see further below)

nadapt <-100 #length of adaptation period.
nburn <- 250 #length of burnin period.
niter <- 100 #number of iterations after burnin. 
nthin <- 1 #frequency of sampling from chain
ncat <- 5 #number of density state categories
nchain <- 4 #number of MCMC chains. 
initsd <- 25 #standard deviation for initial conditions of alpha and beta in MCMC (see alljags section below)
prec <- 0.00001 #precision for prior distributions on alpha and beta
plot.diagnostics <- FALSE #do we make the MCMC diagnostic plots? There are a lot of them, and they aren't here arranged in a particularly convenient way
runMCMC <- TRUE #run the MCMC?
fitML <- FALSE #also do ML estimation as a check?


# Quick check of transition counts
m <- transition.count(mydata,ncat) #summary table of transition counts
print("observed transition counts")
print(m)



############################################################
############################################################







#### break down the dataset into subsets depending on source state ####
sub1 <- mydata[mydata$yt==1,]
sub2 <- mydata[mydata$yt==2,]
sub3 <- mydata[mydata$yt==3,]
sub4 <- mydata[mydata$yt==4,]
sub5 <- mydata[mydata$yt==5,]

#### get explanatory variables in matrix form (only 1 variable (rotation) in this case) ### 
sub1x <- as.matrix(sub1[,3:ncol(sub1)])
sub2x <- as.matrix(sub2[,3:ncol(sub2)])
sub3x <- as.matrix(sub3[,3:ncol(sub3)])
sub4x <- as.matrix(sub4[,3:ncol(sub4)])
sub5x <- as.matrix(sub5[,3:ncol(sub5)])


#### Run the jags model ####
# note that running the content of the curly brackets underneath will likely throw warnings() about incomplete adaptation
if(runMCMC){
  
  # fixed effect baseline-category logit models, conditional on each current state
  jags1.out <- jagsout(sub1x,sub1$ytplus,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec)
  jags2.out <- jagsout(sub2x,sub2$ytplus,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec)
  jags3.out <- jagsout(sub3x,sub3$ytplus,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec)
  jags4.out <- jagsout(sub4x,sub4$ytplus,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec)
  jags5.out <- jagsout(sub5x,sub5$ytplus,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec)
  
  jags1.alpha <- as.mcmc.list(jags1.out$alpha)
  jags2.alpha <- as.mcmc.list(jags2.out$alpha)
  jags3.alpha <- as.mcmc.list(jags3.out$alpha)
  jags4.alpha <- as.mcmc.list(jags4.out$alpha)
  jags5.alpha <- as.mcmc.list(jags5.out$alpha)
  
  jags1.beta <- as.mcmc.list(jags1.out$beta)
  jags2.beta <- as.mcmc.list(jags2.out$beta)
  jags3.beta <- as.mcmc.list(jags3.out$beta)
  jags4.beta <- as.mcmc.list(jags4.out$beta)
  jags5.beta <- as.mcmc.list(jags5.out$beta)
  
  #mcmc diagnostics
  print("jags1.alpha")
  print(gelman.diag(jags1.alpha[,-1,])) #skip alpha_1, fixed at 0
  print(effectiveSize(jags1.alpha))
  
  print("jags1.beta")
  print(gelman.diag(jags1.beta[,c(-1,-6),])) #first and 5th elements fixed at 0
  print(effectiveSize(jags1.beta))
  
  print("jags2.alpha")
  print(gelman.diag(jags2.alpha[,-1,])) #skip alpha_1, fixed at 0
  print(effectiveSize(jags2.alpha))
  
  print("jags2.beta")
  print(gelman.diag(jags2.beta[,c(-1,-6),])) #first and 5th elements fixed at 0
  print(effectiveSize(jags2.beta))
  
  print("jags3.alpha")
  print(gelman.diag(jags3.alpha[,-1,])) #skip alpha_1, fixed at 0
  print(effectiveSize(jags3.alpha))
  
  print("jags3.beta")
  print(gelman.diag(jags3.beta[,c(-1,-6),])) #first and 5th elements fixed at 0
  print(effectiveSize(jags3.beta))
  
  print("jags4.alpha")
  print(gelman.diag(jags4.alpha[,-1,])) #skip alpha_1, fixed at 0
  print(effectiveSize(jags4.alpha))
  
  print("jags4.beta")
  print(gelman.diag(jags4.beta[,c(-1,-6),])) #first and 5th elements fixed at 0
  print(effectiveSize(jags4.beta))
  
  print("jags5.alpha")
  print(gelman.diag(jags5.alpha[,-1,])) #skip alpha_1, fixed at 0
  print(effectiveSize(jags5.alpha))
  
  print("jags5.beta")
  print(gelman.diag(jags5.beta[,c(-1,-6),])) #first and 5th elements fixed at 0
  print(effectiveSize(jags5.beta))
  
  if(plot.diagnostics){
    
    plot(as.mcmc.list(jags1.alpha))
    par(mfrow=c(4,2))
    traceplot(jags1.beta)
    densplot(jags1.beta)
    
    plot(as.mcmc.list(jags2.alpha))
    par(mfrow=c(4,2))
    traceplot(jags2.beta)
    densplot(jags2.beta)
    
    plot(as.mcmc.list(jags3.alpha))
    par(mfrow=c(4,2))
    traceplot(jags3.beta)
    densplot(jags3.beta)
    
    plot(as.mcmc.list(jags4.alpha))
    par(mfrow=c(4,2))
    traceplot(jags4.beta)
    densplot(jags4.beta)
    
    plot(as.mcmc.list(jags5.alpha))
    par(mfrow=c(4,2))
    traceplot(jags5.beta)
    densplot(jags5.beta)
    
  }
  
  #compile MCMC output
  alljags <- list(jags1.out,jags2.out,jags3.out,jags4.out,jags5.out)
  
  
  # save as RDS 
  #saveRDS(alljags, file="alljags_FE_model.rds")
  
}


### Sampling output - Alljags ###


#Alljags is a nested list of parameter estimates, alpha (intercept) & beta (slope).
str(alljags)
# dimensions are [[source state]]$parameter[destination state, iteration, chain].
# so

alljags[[1]]$alpha[1,1,1]
 
# is the estimate for source state 1, to destination state 1, iteration 1 from chain 1 etc...


# construct.p  takes a matrix of sites*explanatory variables (x) to create the transition probability matrix.
x<-as.matrix(mydata[,3])
str(x)


# construct.p creates transition matrices for all mcmc reps and sites. 
P <- construct.p(alljags,x)
#P is an array of transition probability matrices, with dimensions [destination state, source state, site, iteration,chain].
# Thus P[,,2,7,3] is the transition probability matrix for the second site, 7th MCMC iteration, 3rd chain. 




# this is Remi's version of Matt Spencer's much more complex plotP() function, which was designed to deal with continuous explanatory variables
# this creates the transition prob matrices for each category of our explanatory variable 
plotPr1<-function(P,iters,mydata,sourcecat,gosource,fact.r,rotation){
  subdata<-mydata[mydata$yt==sourcecat,] #subset of raw data in which the source state is the one of interest
  Psub<-apply(P[gosource,sourcecat,mydata$yt==sourcecat,iters,],1,mean,na.rm=TRUE)
  unique(Psub[which(subdata[,fact.r]==rotation)])
}


# set up matrices for estimates for all state-to-state transitions 
res.mat<-matrix(0,ncol=5,nrow=5) # no rotation - Second year columns, first year rows. Numbers are density states. 
res.mat1<-matrix(0,ncol=5,nrow=5) # rotation

#  calculating transition matrices for rotated and un rotated fields. 
for(sourcecat in 1:5){
  res.mat[,sourcecat]<-sapply(c(1:5),function(i)plotPr1(P,c(1:niter),mydata,sourcecat,i,3,0))
  res.mat1[,sourcecat]<-sapply(c(1:5),function(i)plotPr1(P,c(1:niter),mydata,sourcecat,i,3,1))
  print(sourcecat)}

# have a look at the matrices....
par(mfrow=c(1,2))
image(t(res.mat1))
image(t(res.mat))




