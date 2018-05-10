#-------------------------------------------- ########## Summary statistics functions ############# ----------------------------------------------#


#-------------------------------------------- #################################################### ----------------------------------------------#

##### Calculates damping ratio for of a matrix #####
damping.ratio <- function(mat){
  
  # mat = transition matrix
  
  ev <- eigen(mat)# Take eigen values from a matrix
  dr <- rle(round(Mod(ev$values), 5))$values # Take dominant and sub-dominant eigen vals
  dr <- dr[1]/dr[2]# Calculate ratio
  return(dr)
}

##### Calculate the stable state distribution of a matrix ####
stable.states.distr <- function(mat){
  
  # mat = transition matrix
  
  ev <- eigen( mat ) # Get eigen values
  
  dist <- ev$vector[,1] / sum(ev$vector[,1]) # Solve system
  
  return(dist) 
}


#### Function to calculate the mean density state & variance ####
mean_ds <- function(ds,states){
  
  
  # ds = counts of how many observations in each density state OR the density state distribution as a proportion
  # states = vector of values for each state 
  
  # get proportions
  fx <- ds/sum(ds)
  
  ms <- fx%*%states # get mean density state
  
  # states <- c(1,2,3,4,5)
  
  vs <- fx%*%((states-c(ms)))^2 # get variance associated with distribution
  
  return(data.frame(ms=ms,vs=vs))
  
}


##### Integrates gamma function #####
gamma.int.fun <- function( breaks, shape, rate ) {
  
  # breaks = the break points that define the numerical abundances contained within each density state
  # shape  = the shape parameter of the gamma distribution to be fitted to the stable density structure
  # rate   = the rate parameter of the gamma distribution to be fitted to the stable density structure
  # both rate and shape are determined using optim below
  
  n <- length(breaks) # Number of density states
  starts <- breaks[1:(n-1)] 
  ends <-  breaks[ 2:n]  
  probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate) # calculates definite integral for a set of shape & rate parameters 
  return( probs )
}

##### Optimizes the gamma function for data and given breaks #####
opt.gamma.fun <- function( ds, breaks, init = c(1, 1) ) {
  
  # --- Function that is optimised, crudely based on least squares --- #
  optFun <- function( pars) {
    
    # takes values for shape & rate parameters of gamma distribution
    # returns sums of squares value  - to be minimised using optim below 
    
    shape <- pars[1]
    rate <-pars[2]
    probs <- gamma.int.fun( breaks, shape, rate)
    ssq <-  ( ds - probs )^2 #  residuals squared 
    
    return( sum( ssq) )
  }
  
  ds <- ds / sum(ds) # proportion of counts in each density category 
  op <- optim( init, optFun) # finds the best fit of shape and rate parameters
  n <- length(breaks) 
  shape <- op$par[1] 
  rate <- op$par[2]
  starts <-  breaks[1:(n-1)] 
  ends <-  breaks[ 2:n] 
  
  probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate) # integral  
  
  
  ## Code to plot 
  
  # x <- seq(min(breaks), max(breaks), by = (max(breaks) - min(breaks)) / 100)
  # 
  # fx <-  pgamma(x[2:length(x) ], shape, rate) - pgamma(x[1:(length(x)-1)], shape, rate)
  # 
  # plot( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, probs, type = "l" , log = "x",
  #       xlab = "State", ylab = "Proportion", ylim = c(0,1))
  # 
  # points( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, ds, col= "red" )
  
  mean <- shape / rate # Average density per sample quadrat?
  var <- shape / rate^2 # Variance
  
  
  return(data.frame(mean=mean,variance=var))
  #return( list(shape = shape, rate = rate, mean = mean, var = var, probs = probs, ds = ds) )
  
}


### All together ###
calc.pop.metrics <- function(field.matrices, breaks){
  
  # Takes:
  # field.matrices - average transition matrices for each field
  # x - matrix of explanatory variables from DS data
  # breaks - break points for numerical abundances in each density category
  
  # Gives:
  # Dataframe of field level damping ratios, average densities and associated variances for each field by rotation 
  
  
  nx <- length(field.matrices) # Number of explanatory vars
  rotation.list <- vector("list", length = nx)  # List for metrics per rotation
  
  try(
    
    for (i in 1:nx){ # split by rotation
      print(i)
      subdat <- field.matrices[[i]][,,,-2] # take rotation i and remove sd matrices
      fn <-  ifelse(is.na(dim(subdat)[3]), 1,dim(subdat)[3])  # get nfields for this rotation
      
      
      pop.metrics.data <- data.frame(matrix(NA,nrow=fn,ncol=4)) # empty dataframe for metrics
      names(pop.metrics.data) <- c("rotation", "DR","meanDens","varDens") # names...
      
      ## Implement getting field ID for each rotational subset...
      #names(pop.metrics.data) <- c("rotation", "field", "DR","meanDens","varDens") # names...
      # fi    <- as.numeric(attr(subdat,which ="dimnames")[[3]]) # vector of field IDs in this subset
      # pop.metrics.data[,"field"]  <- fi # add field IDs to data
      
      pop.metrics.data[,"rotation"] <- names(field.matrices)[i] # add rotation type to data
      
      
      ## Calculate damping ratio for each field ##     
      pop.metrics.data[,"DR"] <- apply(subdat,3,damping.ratio)
      
      ## Calculate stable state distribution for each field ##
      ssd <- apply(subdat,3,stable.states.distr)
      
      ## Calcuate mean density and variance from these stable state distributions ##
      densDat  <- do.call(rbind,apply(ssd,2,opt.gamma.fun, breaks=breaks))
      
      pop.metrics.data[,"meanDens"] <- densDat[,"mean"] # add means
      pop.metrics.data[,"varDens"] <- densDat[,"variance"] # variance
      
      
      rotation.list[[i]] <-  pop.metrics.data # slot into matrix
    }
  )
  
  return(do.call(rbind ,rotation.list))
}



