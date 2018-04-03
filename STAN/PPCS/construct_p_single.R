construct.p.single <- function(random,mydata){
  
  
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
    for(f in 1:length(unique(raneffID))){
      
      eta[raneffID==f,] <-  (as.matrix(x[raneffID==f,]) %*% beta) + g[f]
      
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
