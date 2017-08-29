# code below was part of a side analysis we were doing to determine the best way of analysing the data

# we were looking at three approaches, yours and other multinom models using VGAM and BRMS packages in R

rm(list=ls())


#setwd("C:/Users/Robert/Google Drive/PhD/R files/DS_model_weeds/")
setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/")


library(nnet) 
library(mgcv)
library(coda)
library(MASS)
library(dplyr)
library(rjags)
library(tidyr)
library(ggplot2)




###################### rJags model with test data #############################

###############################################################################



# source needed functions

source("simulate.density.random.R")

source("simulate.densitystructure.R")

#source("chainfuncsmulti.R")

source("construct.P.R")

source("transition.count.R")

source("Psummary.R")

source("plotP.R")

source("chainfuncsmulti_test.R")

source("construct.P.random.R")

source("simulate.density.random.coeff.R")




################# Data setup ########################

# read in
mydata <- read.csv("rotation_2007_2008_2009_2010.csv")

str(mydata)

# select 4 most common rotations
mydata <- na.omit(droplevels(mydata)) %>% 
  filter( Crop1 %in% c("barley","beans","wheat","osr"),
          Crop2 %in% c("barley","beans","wheat","osr"),
          Survey %in% c(2009,2008))



# Check which rotations occur
table(mydata$Crop1, mydata$Crop2)


# get the type of rotations in a single column 
mydata$rt <- interaction(droplevels(mydata$Crop1),droplevels(mydata$Crop2))

# setup model matrix for types of interactions alongside density states and fields
rt.rec <- data.frame(model.matrix(~rt-1, data = mydata))

str(mydata)

# Select columns where rotation occurs
# rt.rec.ind <- sapply(rt.rec, function(col) length(unique(col)))
# rt.rec <- rt.rec[,rt.rec.ind>1]


# select rotations that occured in >=5 fields 
n.field.tab <- as.data.frame.matrix(table(mydata$field,mydata$rt))
n.fields <- sapply(n.field.tab, function(col) length(unique(col)))

rt.rec <- rt.rec[,n.fields>=5]


# create new field index for each survey year
mydata$fieldyear <- interaction(droplevels(mydata$field),droplevels(factor(mydata$Survey)))


mydata <- droplevels(cbind(mydata,rt.rec)) 



  
# setup for jags models 
mydata <- droplevels(mydata) %>% 
  
  
  mutate(yt = Year1 +1 , 
         
         ytplus = Year2 +1,
         
         field= as.integer(fieldyear)) %>%
  
 select (yt, ytplus, field,12:ncol()) 

str(mydata)

head(mydata)

table(mydata$field)

#test numbers: Matt Spencer mentions that are nowhere near large enough but they give reasonable transition estimates, and much larger numbers become unworkable with the number of observations we've got (because of the size of the P array, see further below)
nadapt <-100 #length of adaptation period, meant to be shorter than burnin, and burnin itself should be longer

nburn <- 100#length of burnin period. 

niter <- 100 #number of iterations after burnin. using 250 for spdata1 is too much, the array size becomes unmanageable for building an array 

nthin <- 1 #frequency of sampling from chain

ncat <- 5 #number of abundance categories

nchain <- 4 #number of MCMC chains.

initsd <- 5 #standard deviation for initial conditions of alpha and beta in MCMC

vf.mu <- 4 # shape of vfield gamma prior, also used to set initial conditions for mcmc chains

vf.sd <- 1 # rate of vfield gamma prior, also used to set initial conditions for mcmc chains 

prec <- 0.18 #precision for prior distributions on alpha and beta

nfield <- length(unique(mydata$field))

field.vmax <- 15 #max variance for field effect

plot.diagnostics <- FALSE #do we make the MCMC diagnostic plots? There are a lot of them, and they aren't here arranged in a particularly convenient way

runMCMC <- TRUE #run the MCMC?

fitML <- FALSE #also do ML estimation as a check?


########### Checking with simulated data ##############
# nfield <- 10
# nsim <- 1e2#number of simulated observations: with only 100, we may not do that well...
# ncat <- 5 #number of simulated categories
# pinit <- c(0.2,0.2,0.2,0.2,0.2) #initial state probabilities for simulated data
# alpha <- t(matrix(c(0,  0,  0,  0,  0,
#                   4,4,4,4,4,
#                   3,3,3,3,3,
#                   2,2,2,2,2,
#                   1,1,1,1,1), nrow=ncat,ncol=ncat))
# 
# 
# beta <- array(0.1*1:32,c(5,5,2))
# beta[1,,] <- 0 #identifiability constraint
# 
# vfield.source <- c(0,0,0,0,0) # vfield dependent on source states 1:5 respectively
# 
# 
# mydata <- simulate.random(n = nsim,nfield = nfield ,ncat = ncat,pstart = pinit,alpha = alpha,vfield.source=vfield.source)
# 
# 
# mydata <- simulate.random.coeff(n = nsim,nfield = nfield ,ncat = ncat,pstart = pinit,alpha = alpha,beta=beta,vfield.source=vfield.source)
##################################




# File name for saving model output # 
SaveSamples <-  "alljags_fieldrandom_all_NG.rds"
SaveMods    <- "allmods_fieldrandom_all_NG.rds"




# break down the dataset into subsets depending on source state

sub1 <- mydata[mydata$yt==1,]

sub2 <- mydata[mydata$yt==2,]

sub3 <- mydata[mydata$yt==3,]

sub4 <- mydata[mydata$yt==4,]

sub5 <- mydata[mydata$yt==5,]



# get fields in matrix form (rotation model)

sub1f <- as.matrix(sub1[,3])

sub2f <- as.matrix(sub2[,3])

sub3f <- as.matrix(sub3[,3])

sub4f <- as.matrix(sub4[,3])

sub5f <- as.matrix(sub5[,3])

# get explanatory variables in matrix
sub1x <- as.matrix(sub1[,4:ncol(sub1)])
sub2x <- as.matrix(sub2[,4:ncol(sub2)])
sub3x <- as.matrix(sub3[,4:ncol(sub3)])
sub4x <- as.matrix(sub4[,4:ncol(sub4)])
sub5x <- as.matrix(sub5[,4:ncol(sub5)])




if(runMCMC){
  
  #baseline-category logit models, conditional on each current state
  
  jags1.out <- jagsfieldout.x(field = sub1f,x = sub1x,y = sub1$ytplus,ncat,nchain,nadapt,nburn,niter,nthin,initsd,prec,nfield)
  
  jags2.out <- jagsfieldout.x(field = sub2f,x = sub2x,y = sub2$ytplus,ncat,nchain,nadapt,nburn,niter,nthin,initsd,prec,nfield)
  
  jags3.out <- jagsfieldout.x(field = sub3f,x = sub3x,y = sub3$ytplus,ncat,nchain,nadapt,nburn,niter,nthin,initsd,prec,nfield)
  
  jags4.out <- jagsfieldout.x(field = sub4f,x = sub4x,y = sub4$ytplus,ncat,nchain,nadapt,nburn,niter,nthin,initsd,prec,nfield)
  
  jags5.out <- jagsfieldout.x(field = sub5f,x = sub5x,y = sub5$ytplus,ncat,nchain,nadapt,nburn,niter,nthin,initsd,prec,nfield)
  
  
  

  jags1.alpha <- as.mcmc.list(jags1.out$samples$alpha)

  jags2.alpha <- as.mcmc.list(jags2.out$samples$alpha)

  jags3.alpha <- as.mcmc.list(jags3.out$samples$alpha)

  jags4.alpha <- as.mcmc.list(jags4.out$samples$alpha)

  jags5.alpha <- as.mcmc.list(jags5.out$samples$alpha)

  jags1.vfield <- as.mcmc.list(jags1.out$samples$vfield)

  jags2.vfield <- as.mcmc.list(jags2.out$samples$vfield)

  jags3.vfield <- as.mcmc.list(jags3.out$samples$vfield)

  jags4.vfield <- as.mcmc.list(jags4.out$samples$vfield)

  jags5.vfield <- as.mcmc.list(jags5.out$samples$vfield)



  #mcmc diagnostics

  print("jags1.alpha")

  print(gelman.diag(jags1.alpha[,-1,])) #skip alpha_1, fixed at 0

  print(effectiveSize(jags1.alpha))



  print("jags2.alpha")

  print(gelman.diag(jags2.alpha[,-1,])) #skip alpha_1, fixed at 0

  print(effectiveSize(jags2.alpha))



  print("jags3.alpha")

  print(gelman.diag(jags3.alpha[,-1,])) #skip alpha_1, fixed at 0

  print(effectiveSize(jags3.alpha))



  print("jags4.alpha")

  print(gelman.diag(jags4.alpha[,-1,])) #skip alpha_1, fixed at 0

  print(effectiveSize(jags4.alpha))



  print("jags5.alpha")

  print(gelman.diag(jags5.alpha[,-1,])) #skip alpha_1, fixed at 0

  print(effectiveSize(jags5.alpha))



  print("jags1.vfield")

  print(gelman.diag(jags1.vfield))

  print(effectiveSize(jags1.vfield))



  print("jags2.vfield")

  print(gelman.diag(jags2.vfield))

  print(effectiveSize(jags2.vfield))



  print("jags3.vfield")

  print(gelman.diag(jags3.vfield))

  print(effectiveSize(jags3.vfield))



  print("jags4.vfield")

  print(gelman.diag(jags4.vfield))

  print(effectiveSize(jags4.vfield))



  print("jags5.vfield")

  print(gelman.diag(jags5.vfield))

  print(effectiveSize(jags5.vfield))







    if(plot.diagnostics){



      plot(as.mcmc.list(jags1.alpha))

      plot(as.mcmc.list(jags2.alpha))

      plot(as.mcmc.list(jags3.alpha))

      plot(as.mcmc.list(jags4.alpha))

      plot(as.mcmc.list(jags5.alpha))
      
      plot(as.mcmc.list(jags5.out$samples$beta))

      plot(as.mcmc.list(jags1.vfield))

      plot(as.mcmc.list(jags2.vfield))

      plot(as.mcmc.list(jags3.vfield))

      plot(as.mcmc.list(jags4.vfield))

      plot(as.mcmc.list(jags5.vfield))



    }

  
  
  
  #compile MCMC output
  
  alljags <- list(jags1.out$samples,jags2.out$samples,jags3.out$samples,jags4.out$samples,jags5.out$samples)
  
  allmods <- list(jags1.out$model,jags2.out$model,jags3.out$model,jags4.out$model,jags5.out$model)
  
  #save the MCMC output
  
  
  
  saveRDS(alljags, file="alljags_fieldrandom_all_NG.rds")
  
  
}



f <-  as.matrix(mydata[,3])
x <- as.matrix(mydata[,4:ncol(mydata)])


saveRDS(x, file="jags_x.rds")
saveRDS(f, file="jags_f.rds")



