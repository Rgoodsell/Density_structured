
##### Constructs a single column of the transition matrix #####
construct.pcol.random<-function(alpha,beta,gmma,x,ncat){
  #Calculates one column of the transition probability matrix for a specific mcmc rep and site.
  #inputs
  #alpha is vector of intercepts for this source state (one per destination state, with the first one constrained to be zero)
  #beta is a matrix of slopes for this source state (element i,j is the effect of the jth explanatory variable on the ith destination state)
  #x is a matrix (sites x variables) of explanatory variables
  #ncat is the number of states
  # gma is a matrix of field effects for this source state (element i,j is the random field effect of the ith field on the jth destination state)
  #returns
  #a single column of the transition probability matrix
  
  
  n<-dim(x)[1] # get number of observations
  a1 <-matrix(alpha,nrow=n,ncol=ncat,byrow=TRUE) # get the estimates of alpha from a single site from a single mcmc rep

  phi<-exp(a1 + x%*%t(beta) + gmma) # calculate phi and a column of the transition matrix
  pcol<-phi/rowSums(phi) 
  return(t(pcol))
}

##### Constructs entire transition matrices for each site & mcmc rep #####
construct.p.random<-function(alljags,x,f,nmcmc=NULL){
  #Calculates transition probability matrices for all mcmc reps and sites.
  #input:
  # alljags is a list of MCMC outputs from each source state
  # x is a matrix (sites x variables) of explanatory variables
  # nmcmc is number of mcmc reps to use (uses all if not specified)
  # output
  #P is an array of transition probability matrices, with dimensions destination state, source state, site, MCMC iteration,chain. Thus P[,,2,7,3] is the transition probability matrix for the second site, 7th MCMC iteration, 3rd chain
  
  
## Get the number of...
  nsite<-dim(x)[1] # Observations
  nx <-  dim(x)[2] # Explanatory variables
  
  nchain<-dim(alljags[[1]]$alpha)[3] # Chains
  ncat<-dim(alljags[[1]]$alpha)[1] # Density states
  nfield <- dim(alljags[[1]]$gamma)[1]  # Fields
  
  if(missing(nmcmc) | is.null(nmcmc)){
    nmcmc <- dim(alljags[[1]]$alpha)[2] # MCMC iterations
  }
  
  
  
## P is an array of transition probability matrices, with dimensions destination state, source state, site, MCMC iteration,chain. 
## Thus P[,,2,1,7,3] is the transition probability matrix for the 2nd site, 1st field, 7th MCMC iteration, 3rd chain
   P<-array(dim=c(ncat,ncat,nsite,nmcmc,nchain))
  

  
# This function takes  a really long time to construct these matrices, progress bar just monitors progress...  
#  pb <-   txtProgressBar(min = 0, max = (nchain*niter) ,style = 3 )
#  counter <- 0
  
  
  for(cha in 1:nchain){#this is slow, but don't know how to use an apply function over mcmc lists

    for(r in 1:nmcmc){
 #     counter <-  counter + 1
  #    setTxtProgressBar(pb, counter)
      
      for(i in 1:ncat){
        
        ## Get...
        alpha<-as.mcmc.list(alljags[i][[1]]$alpha) # Alpha estimates for source state i
        a<-alpha[[cha]][r,] # Chain cha and iteration r
        aa<-matrix(a,nrow=1,ncol=ncat) 
        
        beta<-as.mcmc.list(alljags[i][[1]]$beta) # Beta estimates for source state i
        b<-beta[[cha]][r,] # Chain cha, iteration r
        bb<-matrix(b,nrow=ncat,ncol=nx) 
        
        gam<-as.mcmc.list(alljags[i][[1]]$gamma) # Gamma estimates for source state i
        g<-gam[[cha]][r,] # Chain cha and iteration r
        gm <- matrix(g,nrow=nfield,ncol=ncat,byrow = F) # As a matrix corresponding to field ID
        gg <- matrix(nrow=nsite, ncol=ncat)  
        
        for (l in 1:nsite){
          gg[l,] <- gm[f[l],] # As a matrix corresponding to site specific field effect as in x
        }
       
      
        P[,i,,r,cha]<-construct.pcol.random(aa,bb,gg,x,ncat) # Construct....
        
        
     }
    }
   }

 
  
  return(P)
}

##### Calculates field level average matrices & SDs from P array #####
calc.field.matrices <-  function(P,x,f){
  
  ## Takes:
  # P - an array of all transition matrices for each iteration*site*chain combination 
  # x a matrix of explanatory variables by site in DS data
  # a matrix of field IDs by site in DS data
  
  
  nx <- dim(x)[2]  # No. explantatory vars

  rotation.list <- vector("list",length = nx) # empty list for constructs
  names(rotation.list) <- colnames(x) # set names of rotations
  
  
  
  for (rotation in 1:nx){ # For each explanatory variable
    
   # Get... 
    xi <- x[,rotation]        # Sites that underwent rotation x    
    fs <- unique(f[xi==1])    # field IDS of fields in this subset
    fn <- length(fs)          # number of fields in this subset

    # construct array for field level transition matrices to go into..
    # dimensions f.array[,,,1:2]  are the   mean and   standard deviation of all transition matrices for a particular field f.array[,,field,] respectively
    f.array <- array(dim=c(ncat,ncat, fn,2))
    dimnames(f.array) <- list(NULL,NULL,fs,NULL) # Name the field level matrices by field ID
    
      for(field in 1:fn){ # For each field in this set 
      
      xii <- P[,,f==fs[field],,]   # Get all transition matrices for a specific field in the set fs
      f.array[,,field,1] <- apply(xii,1:2,mean) # Calculate the average matrix of these matrices
      f.array[,,field,2] <-  apply(xii,1:2,sd) # SD for transition probabilities

                         }
    
      rotation.list[[rotation]] <-f.array # Place field level matrix array into list
    
                        }
  
  return(rotation.list)
}
