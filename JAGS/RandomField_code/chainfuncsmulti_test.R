jagsout<-function(x,y,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec){

#run the baseline-category logit MCMC for a single source category

#input:

#x is explanatory variables (as a matrix, one row per site, one col per variable)

#y is states next year

#ncat is number of states

#nchain is number of MCMC chains to run

#nadapt is number of adaptation steps prior to burnin

#nburn is number of burnin iterations

#niter is number of sampling iterations

#thin is thinning (e.g. if thin=20, we use every 20th iteration)

#initsd is standard deviation for initial conditions (make this large for overdispersion)

#prec is precision for priors on alpha and beta (e.g. 0.00001)

#output: the object produced by jags.samples from the MCMC run

  



  

  

  

  nx<-dim(x)[2]

 priormvn<-makeprior(ncat,prec) #prior mean and precision for alpha and each col of beta

  jm <- jags.model('marclim.logistic.multi.bug',data = list('x' = x,'y' = y,'N' = dim(x)[1],'K' = ncat,'nx'= nx,a.mu=priormvn$a.mu,a.T=priormvn$a.T),

  inits<-function(){#overdispersed initial conditions: annoyingly, this function can't take any useful arguments

	sd<-initsd

	alpha <- c(NA,rnorm(ncat-1,0,sd))

	beta<-matrix(rnorm(ncat*nx,0,sd),nrow=ncat,ncol=nx)

	beta[1,] <- NA

	return(list(alpha=alpha,beta=beta))

},

  n.chains = nchain,n.adapt = nadapt)



print(list.samplers(jm)) #display the samplers used

  

#this works (no initial condition specified):  

  #jm <- jags.model('marclim.logistic.multi.bug',data = list('x' = x,'y' = y,'N' = dim(x)[1],'K' = ncat,'nx'= nx,a.mu=prioralpha$a.mu,a.T=prioralpha$a.T),

   #n.chains = nchain,n.adapt = nadapt)





  update(jm, nburn) #burnin

  jms <- jags.samples(jm,c('alpha','beta'),n.iter=niter,thin=nthin) #sampling

  return(jms)

}



makeprior<-function(ncat,prec){

  #make a prior mean and covariance matrix for alpha, or for a single column of beta

  #input ncat, the number of categories, and prec, the precision

  #output a list containing a.mu (prior mean vector) and a.T (prior precision matrix) for a multivariate normal prior

  a.mu<-rep(0,ncat-1)

  a.T<-diag(ncat-1)*prec

  return(list(a.mu=a.mu,a.T=a.T))

}



jagsfieldout <- function(field,y,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec, field.vmax, nfield){

#run the baseline-category logit MCMC for a single source category, with random field effect

#input:

#field is a vector of field codes (integers)

#y is states next year

#ncat is number of states

#nchain is number of MCMC chains to run

#nadapt is number of adaptation steps prior to burnin

#nburn is number of burnin iterations

#niter is number of sampling iterations

#thin is thinning (e.g. if thin=20, we use every 20th iteration)

#initsd is standard deviation for initial conditions (make this large for overdispersion)

#prec is precision for priors on alpha and beta (e.g. 0.00001)

#field.vmax is max value for uniform prior on among-field variance

#nfield: number of fields    

#output: the object produced by jags.samples from the MCMC run


  


    priormvn <- makeprior(ncat,prec) #prior mean and precision for alpha

    jm <- jags.model('fieldrandom_all.bug',data = list('field' = as.vector(field),'y' = y,'N' = length(y),'K' = ncat, a.mu=priormvn$a.mu,a.T=priormvn$a.T, 'nfield' = nfield),

                     inits<-function(){#overdispersed initial conditions: annoyingly, this function can't take any useful arguments

                         sd <- initsd

                         alpha <- c(NA,rnorm(ncat-1,0,sd))

                         vfield <- rlnorm(1,4,1)

                         return(list(alpha=alpha, vfield=vfield))

                     },

                     n.chains = nchain,n.adapt = nadapt)



    print(list.samplers(jm)) #display the samplers used

  

    update(jm, nburn) #burnin

    jms <- jags.samples(jm,c('alpha','vfield','gamma'),n.iter=niter,thin=nthin) #sampling


    out.list <- list("model" = jm, "samples" =jms)
    return(out.list)

}




jagsfieldout.x<-function(field,x,y,ncat,nchain,nadapt,nburn,niter,thin,initsd,prec,nfield){
  
  #run the baseline-category logit MCMC for a single source category
  
  #input:
  
  #x is explanatory variables (as a matrix, one row per site, one col per variable)
  
  #y is states next year
  
  #ncat is number of states
  
  #nchain is number of MCMC chains to run
  
  #nadapt is number of adaptation steps prior to burnin
  
  #nburn is number of burnin iterations
  
  #niter is number of sampling iterations
  
  #thin is thinning (e.g. if thin=20, we use every 20th iteration)
  
  #initsd is standard deviation for initial conditions (make this large for overdispersion)
  
  #prec is precision for priors on alpha and beta (e.g. 0.00001)
  
  #output: the object produced by jags.samples from the MCMC run
  
  
  
  
  
  
  
  
  
  
  
  nx<-dim(x)[2]
  
  priormvn<-makeprior(ncat,prec) #prior mean and precision for alpha and each col of beta
  
  jm <- jags.model('fieldrandom_all_coeff.bug',data = list('field' = as.vector(field),'x' = x,'y' = y,'N' = dim(x)[1],'K' = ncat,'nx'= nx,a.mu=priormvn$a.mu,a.T=priormvn$a.T ,'nfield' = nfield, vf.mu=vf.mu, vf.sd=vf.sd),
                   

                   inits<-function(){#overdispersed initial conditions: annoyingly, this function can't take any useful arguments
                     
                     sd<-initsd
                     
                     alpha <- c(NA,rnorm(ncat-1,0,sd))
                     
                     beta<-matrix(rnorm(ncat*nx,0,sd),nrow=ncat,ncol=nx)
                     
                     beta[1,] <- NA
                     
                     vfield <- rlnorm(1,vf.mu, vf.sd)
                     
                     return(list(alpha=alpha,beta=beta,vfield=vfield))
                     
                   },
                   
                   n.chains = nchain,n.adapt = nadapt)
  
  
  
  print(list.samplers(jm)) #display the samplers used
  
  
  
  #this works (no initial condition specified):  
  
  #jm <- jags.model('marclim.logistic.multi.bug',data = list('x' = x,'y' = y,'N' = dim(x)[1],'K' = ncat,'nx'= nx,a.mu=prioralpha$a.mu,a.T=prioralpha$a.T),
  
  #n.chains = nchain,n.adapt = nadapt)
  
  
  
  
  
  update(jm, nburn) #burnin
  
  jms <- jags.samples(jm,c('alpha','beta','vfield','gamma'),n.iter=niter,thin=nthin) #sampling
  
  
  out.list <- list("model" = jm, "samples" =jms)
  
  
  
  return(out.list)
  
}



makeprior<-function(ncat,prec){
  
  #make a prior mean and covariance matrix for alpha, or for a single column of beta
  
  #input ncat, the number of categories, and prec, the precision
  
  #output a list containing a.mu (prior mean vector) and a.T (prior precision matrix) for a multivariate normal prior
  
  a.mu<-rep(0,ncat-1)
  
  a.T<-diag(ncat-1)*prec
  
  return(list(a.mu=a.mu,a.T=a.T))
  
}


