construct.pcol<-function(alpha,beta,x,ncat){
#Calculates one column of the transition probability matrix for a specific mcmc rep and site.
#inputs
#alpha is vector of intercepts for this source state (one per destination state, with the first one constrained to be zero)
#beta is a matrix of slopes for this source state (element i,j is the effect of the jth explanatory variable on the ith destination state)
#x is a matrix (sites x variables) of explanatory variables
#ncat is the number of states
#returns
#a single column of the transition probability matrix
	n<-dim(x)[1]
#	phi<-array(0,c(n,ncat))
        a<-matrix(rep(alpha,n),nrow=n,ncol=ncat,byrow=TRUE)
        phi<-exp(a+x%*%t(beta))
#	for(i in 1:ncat){
#		phi[,i]<-exp(alpha[i]+x%*%beta[i,])
#	}

 	pcol<-phi/rowSums(phi) #we will get NaN entries if any infinite phi
        #pcol<-checkinf.pcol(phi) #if only one infinite phi: but decided this is not the right thing to do
	return(t(pcol))
}

#calculate p_ijm from phi. If only one infinite phi in a row, can we set the corresponding transition probability to 1? Concluded that we can't do this, because each transition probability is a function of more than one explanatory variable, and thus in general its limit as all the explanatory variables go to infinity does not exist. Therefore don't use this function.
checkinf.pcol<-function(phi){
  pcol<-phi/rowSums(phi)
  winf<-!is.finite(phi)
  ninf<-rowSums(winf)
  winf[ninf>1,]<-FALSE 
  pcol[winf]<-1
  return(pcol)
}

construct.p<-function(alljags,x,nmcmc=NULL){
#Calculates transition probability matrices for all mcmc reps and sites.
#input
#alljags is a list of MCMC outputs from each source state
#x is a matrix (sites x variables) of explanatory variables
#nmcmc is number of mcmc reps to use (uses all if not specified)
#output
#P is an array of transition probability matrices, with dimensions destination state, source state, site, MCMC iteration,chain. Thus P[,,2,7,3] is the transition probability matrix for the second site, 7th MCMC iteration, 3rd chain

 	nsite<-dim(x)[1]
        nchain<-length(as.mcmc.list(alljags[1][[1]]$alpha))
  	ncat<-dim(as.mcmc.list(alljags[1][[1]]$alpha)[[1]])[2]
	if(missing(nmcmc) | is.null(nmcmc)){
		nmcmc <- dim(as.mcmc.list(alljags[1][[1]]$alpha)[[1]])[1]
		}
 	nx<-dim(x)[2]
	P<-array(dim=c(ncat,ncat,nsite,nmcmc,nchain))
        for(cha in 1:nchain){#this is slow, but don't know how to use an apply function over mcmc lists
          for(r in 1:nmcmc){
            for(i in 1:ncat){
              alpha<-as.mcmc.list(alljags[i][[1]]$alpha)
              a<-alpha[r,][[cha]]
              aa<-matrix(a,nrow=1,ncol=ncat)
              beta<-as.mcmc.list(alljags[i][[1]]$beta)
              b<-beta[r,][[cha]]
              bb<-matrix(b,nrow=ncat,ncol=nx)
              P[,i,,r,cha]<-construct.pcol(aa,bb,x,ncat)
            }
          }
        }
	return(P)
}


meanp<-function(alljags,x){
#Calculates transition probability matrices for each site, using posterior mean parameter estimates.
#input
#alljags is a list of MCMC outputs from each source state
#x is a matrix (sites x variables) of explanatory variables
#output
#Pm is an array of transition probability matrices, with dimensions destination state, source state, site. Thus P[,,2] is the transition probability matrix for the second site.
#NOTE: this isn't a very good point estimate. Use it only for rough things (written only to produce a quick figure for BES report). Much better to get the posterior mean of the summary statistics that are nonlinear functions of transition probabilities, that in turn are nonlinear functions of the underlying parameters. Also not checked very carefully!

 	nsite<-dim(x)[1]
  	ncat<-dim(as.mcmc.list(alljags[1][[1]]$alpha)[[1]])[2]
 	nx<-dim(x)[2]
	meanP<-array(dim=c(ncat,ncat,nsite))
        for(i in 1:ncat){
          alpha<-as.mcmc.list(alljags[i][[1]]$alpha)
          a <- summary(alpha)
          aa<-matrix(a$statistics[,1],nrow=1,ncol=ncat)
          beta<-as.mcmc.list(alljags[i][[1]]$beta)
          b <- summary(beta)
          bb<-matrix(b$statistics[,1],nrow=ncat,ncol=nx)
          meanP[,i,]<-construct.pcol(aa,bb,x,ncat)
        }
        return(meanP)
}


