
###### Predict from ordered logit model ########

rm(list=ls())

library(ggplot2)
library(rstan)
library(dplyr)
library(reshape2)
library(gridExtra)


###### Construct mean thetas from mean parameter estimates #####
## Quick and dirty, should only be used for a quick check ######

construct.theta.mean <- function(betas,cutpoints,x,N){

  
## Takes:
  # betas - column vector of mean beta parameters
  # cutpoints - vector of  mean cutpoints
  # x - matrix of source states
  # N - number of observations
  
## Returns:
  # Theta : matrix of transition probabilities at site n
  
# set up vectors for storing outputs  
eta <- matrix(0,nrow=N)  # linear predictor for each site observation
theta <- matrix(0,nrow=N,ncol=K) # N * K matrix of probabilities of transition from states [,1:K] from site [n,]


eta = as.matrix(x) %*% as.matrix(beta) #  construct linear predictor 

theta[,1] = 1- plogis(eta[,1]-c[1]) # construct theta 1

for(k in 2:(K-1)){
  theta[,k] = plogis(eta - c[k-1]) - plogis(eta-c[k]) # theta 2,3 & 4
}

theta[,K] = plogis(eta - c[K-1]) # theta 5

return(theta)

}

################################################################



######### Construct theta from fixed effect model ##############
### Constructs theta for all parameter estimates after warmup ##

ordered.construct.p <- function(fixed,mydata){

  
  ## Takes:
  # fixed - fixed effect orderedLogit  stan model object
  # mydata - data used to run the model
  
  ## Returns:
  # P : matrix of transition probabilities to state [,1:K,] at site[n,,] for iteration [,,i]


# Extract estimates from stanfit object
stanfit<- rstan::extract(fixed)  

# Get test numbers  
N <- nrow(mydata) # Number of sites/observations
K <- length(unique(mydata$yt)) # Number of categories
x <- data.frame(model.matrix( ~yt -1, data = mydata)) # Matrix of source states at site n
niter <-  nrow(stanfit$beta) # Number of iterations

# Set up output objects
eta <- matrix(0,nrow=N)  # linear predictor for each site observation
theta <- matrix(0,nrow=N,ncol=K) # N * K matrix of probabilities of transition from states [,1:K] from site [n,]
P <-  array(0,dim=c(N,K,niter)) # N * K matrix of probabilities of transition from states [,1:K,] from site [n,,] for each iteration across all chains [,,i]



###  Construct probabilites for all sites from parameter estimates from iteration i ###
for(i in 1:(niter)){

# Get betas and cutpoints for iteration i     
beta <- stanfit$beta[i,]
c <- stanfit$c[i,]

eta = as.matrix(x) %*% as.matrix(beta) #  construct linear predictor 

# Construct theta
theta[,1] = 1- plogis(eta[,1]-c[1]) # construct theta 1
for(k in 2:(K-1)){
  theta[,k] = plogis(eta - c[k-1]) - plogis(eta-c[k]) # theta 2,3 & 4
}
theta[,K] = plogis(eta - c[K-1]) # theta 5

# return probabilties from iteration i
P[,,i] <- theta
}

# Return all transition probability estimates across all parameter estimates. 
return(P)

}

################################################################


#### Construct P for ordered logit model with random effect #### 
### Constructs theta for all parameter estimates after warmup ##


ordered.construct.p.random <- function(random,mydata){
  
  
  ## Takes:
  # fixed - fixed effect orderedLogit  stan model object
  # mydata - data used to run the model
  
  ## Returns:
  # P : matrix of transition probabilities to state [,1:K,] at site[n,,] for iteration [,,i]
  
  # Extract estimates from stanfit object
  stanfit<- rstan::extract(random)
  
  # Get test numbers  
  N <- nrow(mydata) # Number of sites/observations
  K <- length(unique(mydata$yt)) # Number of categories
  x <- data.frame(model.matrix( ~yt -1, data = mydata)) # Matrix of source states at site n
  niter <-  nrow(stanfit$beta) # Number of iterations
  raneffID <- as.matrix(as.numeric(mydata[,3])) # Index of random group levels by site/observation
  
  
  # Set up output objects
  eta <- matrix(0,nrow=N)  # linear predictor for each site observation
  theta <- matrix(0,nrow=N,ncol=K) # N * K matrix of probabilities of transition from states [,1:K] from site [n,]
  P <-  array(0,dim=c(N,K,niter)) # N * K matrix of probabilities of transition from states [,1:K,] from site [n,,] for each iteration across all chains [,,i]

  
  
  # Construct probabilites for all sites from parameter estimates from iteration i
  for(i in 1:(niter)){
    
    # Get betas, gammas (random effects) and cutpoints for iteration i     
    beta <- stanfit$beta[i,]
    c <- stanfit$c[i,]
    g <- stanfit$gamma[i,,]
    
    
   
    
    #  construct linear predictor 
    for(f in 1:length(unique(raneffID))){
    
    # Add beta & gamma before matrix multiplication    
    beta.gamma <-  as.matrix(beta) + as.matrix(g[,f]) 
    
    # Construct eta for each set of sites/observations in field f  
    eta[raneffID==f,]= as.matrix(x[raneffID==f,]) %*% beta.gamma
    
    }
    
   
    # Construct theta
    theta[,1] = 1- plogis(eta-c[1]) # construct theta 1
    for(k in 2:(K-1)){
      theta[,k] = plogis(eta - c[k-1]) - plogis(eta-c[k]) # theta 2,3 & 4
    }
    theta[,K] = plogis(eta - c[K-1]) # theta 5
    
    # return probabilties from iteration i
    P[,,i] <- theta
  }
  
  # Return all transition probability estimates across all parameter estimates. 
  return(P)
  
}

################################################################

############### Construct P single gamma ##############
## Constructs thetas for an orderedLogit model with only a single random intercept (i.e. across all destiation states)
# This model performs significantly worse than a multi-state random intercept and this code is only included  only for reference  


ordered.construct.p.sg <- function(random,mydata){
  
  
  ## Takes:
  # fixed - fixed effect orderedLogit  stan model object
  # mydata - data used to run the model
  
  ## Returns:
  # P : matrix of transition probabilities to state [,1:K,] at site[n,,] for iteration [,,i]
  
  # Extract estimates from stanfit object
  stanfit<- rstan::extract(random)
  
  # Get test numbers  
  N <- nrow(mydata) # Number of sites/observations
  K <- length(unique(mydata$yt)) # Number of categories
  x <- data.frame(model.matrix( ~yt -1, data = mydata)) # Matrix of source states at site n
  niter <-  nrow(stanfit$beta) # Number of iterations
  raneffID <- as.matrix(as.numeric(mydata[,3])) # Index of random group levels by site/observation
  
  
  # Set up output objects
  eta <- matrix(0,nrow=N)  # linear predictor for each site observation
  theta <- matrix(0,nrow=N,ncol=K) # N * K matrix of probabilities of transition from states [,1:K] from site [n,]
  P <-  array(0,dim=c(N,K,niter)) # N * K matrix of probabilities of transition from states [,1:K,] from site [n,,] for each iteration across all chains [,,i]
  
  
  
  # Construct probabilites for all sites from parameter estimates from iteration i
  for(i in 1:(niter)){
    
    # Get betas, gammas (random effects) and cutpoints for iteration i     
    beta <- stanfit$beta[i,]
    c <- stanfit$c[i,]
    g <- stanfit$gamma[i,]
    
    
    
    
    #  construct linear predictor 
    for(f in 1:length(unique(field))){
      
      eta[field==f,]= (as.matrix(x[field==f,]) %*% beta) + g[f]
      
    }
    
    
    
    # Construct theta
    theta[,1] = 1- plogis(eta-c[1]) # construct theta 1
    for(k in 2:(K-1)){
      theta[,k] = plogis(eta - c[k-1]) - plogis(eta-c[k]) # theta 2,3 & 4
    }
    theta[,K] = plogis(eta - c[K-1]) # theta 5
    
    # return probabilties from iteration i
    P[,,i] <- theta
  }
  
  # Return all transition probability estimates across all parameter estimates. 
  return(P)
  
  
}

################################################################


######### Predict density state outcomes ##########################
# Predicts density states next year from source state #############
# Requires average/median  thetas on a site by site basis 

predict.ordered <- function(avgP,mydata){
  
  ## Takes:
  # avgP: Mean/median measure of theta for each site/observation
  # mydata: data used to run the model
  
  ## Returns:
  # Dataframe containing:
  # yt - observed source states
  # ytpreds - predicted destination states
  # ytplus - observed destination states
  # ranIndex - index of random intercept groups if present
  
  
  # Set up predictions depending on whether we used a random intercept.
  if(dim(mydata)[2]>2)
  {
    preds <- matrix(0,nrow=N,ncol=4)
    preds[,4] <- mydata[,3]
    
  } else {
    preds <- matrix(0,nrow=N,ncol=3)
  }
  
  preds[,1] <- mydata$yt
  preds[,2] <- mydata$ytplus
  
  
  ## Sample  from cat dist for observation
  for (n in 1:N){
    preds[n,3] <- sample(1:5,1,prob=avgP[n,]) 
  }
  
  ## Return:
  
  if(dim(mydata)[2]>2){
    full <- data.frame(yt=preds[,1],ytplus=preds[,2],ytpred = preds[,3],ranIndex=preds[,4])
  } else{
    full <-  data.frame(yt=preds[,1],ytplus=preds[,2],ytpred = preds[,3])
  }
  
  return(full)
  
  
  
}


################################################################



############### "Confusion" transition matrices #################
# Plots a matrix of predicted/observed transitions between states


confusion.matrix <- function(preds,mydata){
  
  
  ## Takes:
  # preds: observed source states [,1] and predicted  destination states [,3] 
  # mydata: data used to run the model
  
  ## Returns:
  # ggplot transition matrix containing a ratio of Predicted / observed transitions between states
  
  
  # Make transition matrices for observed and prediction matrices
  predsMat <- data.matrix(table(preds[,1],preds[,3]))
  obsMat <- data.matrix(table(mydata$yt,mydata$ytplus))
  
  # Calculate the ratio of predicted to observed transitions
  accMat <- round((predsMat/obsMat),3)
  
  # Plot it
  accTrans  <-  ggplot(melt(accMat), aes(factor(Var2), factor(Var1))) + geom_tile(aes(fill = value),colour = "white") + 
    geom_text(aes(label=value, fontface="bold"), size=6)+
    scale_fill_gradient(low = "White",high = "Red") +
    scale_x_discrete(name="State in Year 1",
                     breaks=c("1", "2", "3","4","5"),
                     labels=c("A", "L", "M","H","VH"))+
    scale_y_discrete(name="State in Year 2",
                     breaks=c("1", "2", "3","4","5"),
                     labels=c("A", "L", "M","H","VH"))+
    theme(axis.title.x = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x = element_text(size=17),
          axis.text.y = element_text(size=17))+ guides(fill=FALSE)
  
  return(accTrans)
  
}

################################################################



########## Plot distribution for each random effect grouping ###########
# Plots observed vs predicted field/farm level DS state distributions #
# Primarily used to assess model fit 

plot.raneff.matrices <- function(P, mydata){
  
  
  ## Takes:
  # fixed - fixed effect orderedLogit  stan model object
  # mydata - data used to run the model
  
  ## Returns:
  # Plots of observed distributions per group vs median predicted distributions and 5/95% quantile interval
  
  
  # Set up test numbers
  K <- length(unique(mydata[,1])) # Number of categories
  Nraneff <- length(unique(mydata[,3])) # Number of random intercepts
  raneffID <- as.matrix(as.numeric(mydata[,3])) # Random intercept index by site/observation
  raneff.distribs <- array(0,dim=c(K,4,Nraneff)) # Array for calculated distribution
  
  
  # Subset by field ID #
  for(f in 1:Nraneff){
    
    
    fsub_preds <-  P[raneffID==f,,] # Select theta estimates by random effect index
    raneff.distribs[,1,f] <- apply(fsub_preds,2,median) # Calculate median for each destination state for all iterations
    raneff.distribs[,2,f] <- apply(fsub_preds,2, quantile,probs=0.05) # 5% quantile
    raneff.distribs[,3,f] <- apply(fsub_preds,2, quantile,probs=0.95) # 95% quantile
    
    
    # Note: model constraints (i.e. K1 = 0 constraint for identifiability) means that fields/farms
    # with high levels of K1 source states have extremely small density intervals.
    # This is likely due to the abscence of random intercept variation introduced from other 
    # source states, as variation from source K1 = 0. 
    # This needs thorough checking. 
    
    

    fsub_obs <- mydata[raneffID==f,]
    obsMat <- data.matrix(table(fsub_obs[,2]))
    obs <- obsMat/sum(obsMat)
  
    raneff.distribs[1:length(obs),4,f] <- obs
    
    plot(raneff.distribs[,1,f],type="b",col="black",bg="red",pch=21,ylim=c(0,1))
    points(raneff.distribs[,4,f],type="b")
    arrows(x0 = c(1,2,3,4,5),y0 = raneff.distribs[,2,f],
            x1 =c(1,2,3,4,5), y1 = raneff.distribs[,3,f],code = 3,angle=90, length = 0.1)

  }
  
  
}

################################################################

