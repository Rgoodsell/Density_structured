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
  xf <- matrix(0,nrow=nfield,ncol=2)
  xf[1:floor(median(1:nfield)),1] <- 1
  xf[floor(median(1:nfield+1)):nfield,2] <- 1
  
  
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
          phi[i,j,f] <- exp(alpha[i,j] + xf[f,]%*%beta[i,j,] + gamma[i,j,f])
          
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
            ytplus[l,3] <- yt[l,3] # rotational covariates
            ytplus[l,4] <- yt[l,4]
          }
          
        }
      }}
    return(data.frame(yt=yt[,1],ytplus=ytplus[,1],field=yt[,2],rt1=ytplus[,3],rt2=ytplus[,4]))
  }
  
 #  #
 #  #
 #  ##### Test Numbers ######
 #  nfield <- 10 #number of fields
 #  nsim <- 1e3 #number of simulated observations
 #  ncat <- 5 #number of simulated categories
 #  pstart<- c(0.5,0.4,0.3,0.2,0.1) #initial state probabilities for simulated data
 # 
 #  ###### alpha #####
 #  alpha <- t(matrix(c(0,  0,  0,  0,  0,
 #                    0.4,0.4,0.4,0.4,0.4,
 #                    0.3,0.3,0.3,0.3,0.3,
 #                    0.2,0.2,0.2,0.2,0.2,
 #                    0.1,0.1,0.1,0.1,0.1), nrow=ncat,ncol=ncat))
 # 
 #  alpha[1,] <- 0 # set constraint
 # 
 # 
 #  yt <- matrix(rep(0,n), ncol=4, nrow=n)
 #  yt[,1] <- sample(c(1:ncat),replace=T,prob=pstart,size=n)
 #  yt[,2] <- sample(c(1:nfield), replace=T,size=n)
 # 
 #  # assign rotational covariates (either 1 or 0) for field, making sure they're consistent (i.e. fields should only have one type of rotation)
 #  # fields 1:5 have rotation type 1, fields 6:10 have rotation type 2
 #  
 #  # assigns rotational covariates depening on whether field ID is odd or even ¯\_(ツ)_/¯
 # 
 #  for (i in 1:n){
 #  if(yt[i,2] %%  2 == 0)
 #    {yt[i,3] <- 1}
 #    else(yt[i,4] <- 1)}
 # 
 # 
 #  x <- cbind(yt[,3],yt[,4])
 # 
 #  # matrix of explanatory variables by field
 #  xf <- matrix(0,nrow=nfield,ncol=nx)
 #  xf[1:5,1] <- 1
 #  xf[5:10,2] <- 1
 #  xf
 # 
 #  ytplus <- matrix(rep(0,n),ncol=4,nrow=n)
 # 
 # 
 #  # here x needs 2 dimensions (per explanatory) instead of n
 #  # achieve this by creating an extra column which indexes field
 #  # need to make sure this matches up with field in yt and yt plus
 #  # within a field should have the same rotation covariate.
 # 
 # 
 # 
 #  beta <- array(0.1*1:32,c(5,5,2))
 #  beta[1,,] <- 0 #identifiability constraint
 # 
 #  ##### Between field variance and gamma ######
 # 
 #  nsim <- 1000
 #  nfield <- 500
 #  nx <- 2
 #  vfield.source <- (1:5)
 #  ncat <- 5
 # 
 #  vfield.source <- c(1,2,3,4,5) #
 # 
 # 
 #  mydata <- simulate.random.coeff(n = nsim,nfield = nfield ,ncat = ncat,pstart=pinit,alpha = alpha,beta=beta,vfield.source=vfield.source)
 # sum(mydata$rt1)
 # sum(mydata$rt2)
 # 
 # 
