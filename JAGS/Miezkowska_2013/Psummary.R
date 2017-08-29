Psummary<-function(P,iters,spdata){
  #calculate a summary of expected transition counts from MCMC output
  #P is an array of transition probability matrices, with dimensions destination state, source state, site, MCMC iteration,chain. Thus P[,,2,7,3] is the transition probability matrix for the second site, 7th MCMC iteration, 3rd chain
  #iters is the iterations we want to average over (e.g.1:100)
  #spdata is the raw data, including column yt, which indexes current state
  #output
  #Ps is a matrix of predicted counts of transitions, averaged across sites, iterations, and chains (removing any NA values)
  #nNA is a matrix of counts of NAs that were removed
  
  ncat<-dim(P)[1]
  Ps<-matrix(nrow=ncat,ncol=ncat)
  nNA<-Ps
  for(i in 1:ncat){
    Ps[,i]<-apply(P[,,spdata$yt==i,iters,],c(1,2),mean,na.rm=TRUE)[,i]*sum(spdata$yt==i)
    nNA[,i]<-apply(P[,,spdata$yt==i,iters,],c(1,2),function(x){sum(is.na(x))})[,i]
  }
  return(list(Ps=Ps,nNA=nNA))
}
