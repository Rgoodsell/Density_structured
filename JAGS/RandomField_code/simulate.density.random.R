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
  
  
  
  
  
  
  
  # #
  # #
  # ##### Test Numbers ######
  # nfield <- 10000 #number of fields
  # nsim <- 1e5 #number of simulated observations
  # ncat <- 5 #number of simulated categories
  # pstart<- c(0.5,0.4,0.3,0.2,0.1) #initial state probabilities for simulated data
  # 
  # ###### alpha #####
  # alpha <- t(matrix(c(0,  0,  0,  0,  0,
  #                   0.4,0.4,0.4,0.4,0.4,
  #                   0.3,0.3,0.3,0.3,0.3,
  #                   0.2,0.2,0.2,0.2,0.2,
  #                   0.1,0.1,0.1,0.1,0.1), nrow=ncat,ncol=ncat))
  # 
  # alpha[1,] <- 0 # set constraint
  # 
  # ##### Between field variance and gamma ######
  # 
  # nfield <- 10
  # vfield.source <- (1:5)
  # ncat=5
  # 
  # gamma=list()
  # gamma<-  replicate(nfield ,diag(5),simplify = F)
  # 
  # 
  # 
  # 
  # vfield.source <- c(1,2,3,4,5) #
  # 
  # 
  # 
  # 
  # data <- simulate.random (nsim,ncat,nfield,pstart,alpha,vfield.source)
  # 
  # str(data)
  # 
  # 
  # 
  # 