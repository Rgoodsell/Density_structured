##### Calculates damping ratio for of a matrix #####
damping.ratio <- function(mat){
  
  
  ev <- eigen(mat)# Take eigen values from a matrix
  dr <- rle(round(Mod(ev$values), 5))$values # Take dominant and sub-dominant eigen vals
  dr <- dr[1]/dr[2]# Calculate ratio
  dr
}

##### Calculate the stable state distribution of a matrix ####
stable.states.distr <- function(matrix){
  
  ev <- eigen( matrix ) # Get eigen values
  
  dist <- ev$vector[,1] / sum(ev$vector[,1]) # Solve system
  
  return(dist) 
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
  
  
  return(data.frame(mean=mean,var=var))
  #return( list(shape = shape, rate = rate, mean = mean, var = var, probs = probs, ds = ds) )
  
}



### All together ###
calc.pop.metrics <- function(field.matrices,x, breaks){
  
  # Takes:
  # field.matrices - average transition matrices for each field
  # x - matrix of explanatory variables from DS data
  # breaks - break points for numerical abundances in each density category
  
  # Gives:
  # Dataframe of field level damping ratios, average densities and associated variances for each field by rotation 
  

nx <- dim(x)[2] # Number of explanatory vars
rotation.list <- vector("list", length = nx)  # List for metrics per rotation
  


  for (i in 1:nx){ # split by rotation

    subdat <- field.matrices[[i]][,,,-2] # take rotation i and remove sd matrices
    fn <- dim(subdat)[3]                 # number of fields in this subset
    
    pop.metrics.data <- data.frame(matrix(NA,nrow=fn,ncol=5)) # empty dataframe for metrics
    names(pop.metrics.data) <- c("rotation", "field", "DR","meanDens","varDens") # names...
    fi    <- as.numeric(attr(subdat,which ="dimnames")[[3]]) # vector of field IDs in this subset
    pop.metrics.data[,"field"]  <- fi # add field IDs to data
    pop.metrics.data[,"rotation"] <- names(field.matrices)[i] # add rotation type to data


    ## Calculate damping ratio for each field ##     
    pop.metrics.data[,"DR"] <- apply(subdat,3,damping.ratio)

    ## Calculate stable state distribution for each field ##
    ssd <- apply(subdat,3,stable.states.distr)
    
    ## Calcuate mean density and variance from these stable state distributions ##
    densDat  <- do.call(rbind,apply(ssd,2,opt.gamma.fun, breaks=breaks))
    
    pop.metrics.data[,"meanDens"] <- densDat[,"mean"] # add means
    pop.metrics.data[,"varDens"] <- densDat[,"var"] # variance


    rotation.list[[i]] <-  pop.metrics.data # slot into matrix
                  }


return(do.call(rbind ,rotation.list))
}


# breaks <- c(0,1,160,450,1450,2000)
# pop.dat <-calc.pop.metrics(field.matrices,nfield, x, breaks)


