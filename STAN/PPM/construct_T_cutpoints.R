########## Construct T. for cutpoints model #########

# Function to split stanfit into a list by iteration - allows vectorisation
split_stanfit_cutpoints <- function(stanfit,pars){
  
  
  #stanfit = stanfit model object
  #pars = vector of parameters needed to construct matrices (beta & cutpoints)
  
  # extract stanfit
  mod_pars <- rstan::extract(stanfit, pars=pars)
  
  # Empty list for output 
  out_list <- list()
  
  # Loop over iteration indices and store as list of parameter values 
  for(i in 1:nrow(mod_pars[[1]])){
    # ifcond in anon fun below to account for array dim var bewtween params 
    b <-  mod_pars[[1]][i,]
    c <- mod_pars[[2]][i,,]
    
    out_list[[i]] <- list(b,c)  
    # Name list
    names(out_list[[i]]) <- pars
    
  }
  # Return nested lists
  return(out_list)
}


# Function to calculate all field-level transition matrices from a single iteration #
construct_T_cutpoints <- function(input.list,nfield, K){
  
  
  #input.list = list of parameters split by iteration
  #nfield = number of group level effects
  #K = number of categories
  
  # Set up output objects
  eta <- array(0,dim=c(K,nfield))  # source state specific linear predictor for each field
  theta <- array(0,dim=c(nfield,K,K)) # Probabilities, row = field, column = densitination, dim3 = source
  
  beta <-  matrix(input.list$beta, ncol=K,nrow=nfield,byrow=T)  
  c <- input.list$c
  
  
  # Construct a liear predictor for that represents all source states
  # In this case just the estimate of beta.. 
  eta <- beta 
  
  # For source states 1:K
  for(s in 1:K){
    # construct transition probability from source state s, to all densitation states, vectorised over field. 
    theta[,1,s] = 1- plogis(eta[,s]-c[,1]) # Destination state 1
    
    
    for(k in 2:(K-1)){
      theta[,k,s] = plogis(eta[,s] - c[,k-1]) - plogis(eta[,s]-c[,k]) # DS 2,3 & 4
    }
    
    theta[,K,s] = plogis(eta[,s] - c[,K-1]) # DS 5
  }
  
  
  # Transform to transition matrices
  T.  <- lapply(split(theta,seq(nrow(theta))),
                function(x) matrix(x,nrow=5,ncol=5))
  
  
  # Return all transition probability matrices for all fields
  return(T.)
  
  
}


# Gets an array of average T. matrices and there SDS for each field
average_matricies_cutpoints <- function(split,nfield,niter,K){
  
  
  # Empty array for matrices
  mean_matrices <- array(0,dim=c(K,K,nfield,2))
  
  # Get field level matrices by mcmc iteration in an array
  iter_list <- lapply(split,construct_T_cutpoints,nfield=nfield,K=K) %>% unlist(.) %>% 
    array(., dim=c(K, K, nfield,niter))
  
  
  # Average over field & take SDs
  mean_matrices[,,,1] <- apply(iter_list,c(1,2,3),mean)
  mean_matrices[,,,2] <- apply(iter_list,c(1,2,3),sd)
  
  
  
  return(mean_matrices)
}



