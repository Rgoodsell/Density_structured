simulate.densitystructure<-function(n,ncat,pstart,alpha,beta){
  #simulate a data set with state yt sampled from a multinomial, and state yt+1 sampled from a baseline-category logit model conditional on yt, with two explanatory variables (whose values are generated from standard normal distributions)
  #n is number of pairs of observations
  #ncat is number of categories
  #pstart is a vector of ncat elements: the parameters of the multinomial distribution for state at time t
  #alpha is an ncat x ncat array, where element alpha_ij is the intercept for category i next year, conditional on category j this year. Identifiability constraint: alpha_0j should equal zero for all j
  #beta is an ncat x ncat x 2 array, where element beta_ijk is the effect of explanatory variable k on the linear predictor for category i next year, conditional on category j this year. Identifiability constraint: beta_0jk should equal 0 for all j and k
  #
  #returns a data frame containing yt,ytplus (simulated density categories at times t,t+1),scalefetch,scalesst (the simulated explanatory variables, sampled from N(0,1),fetch,sst (identical to scalefetch and scalesst: in here only because in the real data, the unscaled explanatory variables don't have standard deviation 1))
  #yt is state this year (sampled from multinomial with probabilities pstart)
  #ytplus is state next year (generated from baseline-category logit model conditional on yt)
  
  yt <- sample(c(1:ncat),replace=T,prob=pstart,size=n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x <- cbind(x1,x2)
  ytplus <- rep(0,n)
  for(j in 1:ncat){#source state
    phi <- array(0,c(n,ncat))
    for(i in 1:ncat){#destination state
      #construct linear predictor for the jth source state and ith destination state
      phi[,i] <- exp(alpha[i,j]+x%*%beta[i,j,])
    }
    p <- phi/rowSums(phi)
    for(l in 1:n){
      if(yt[l]==j){
        ytplus[l] <- sample(c(1:ncat),prob=p[l,],size=1)#sample from categorical distribution with probabilities from baseline-category logit model conditional on current state
      }
    }
  }
  return(data.frame(yt=yt,ytplus=ytplus,scalefetch=x1,scalesst=x2,fetch=x1,sst=x2))
}
