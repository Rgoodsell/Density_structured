#functions for analysis of transition probability matrices e.g. stationary distributions, ergodicity

stat.dist<-function(P){
	# Stationary distribution of a set of transition probability matrices for many sites and mcmc reps.
	# P is an ncat x ncat x nsite x nmcmc x nchain array, where ncat is the no. of categories in the model, nsite is the no. of sites, and nmcmc is the number of mcmc reps and nchain is the number of chains
	# The output is w, an ncat x nsite x nmcmc x nchain array whose columns are stationary distributions for given sites, mcmc reps and chains, and dr, nsite x nmcmc x nchain array whose elements are damping ratios for given sites, mcmc reps, and chains
	ncat<-dim(P)[1]
	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
 	w<-array(dim=c(ncat,nsite,nmcmc,nchain))
        dr<-array(dim=c(nsite,nmcmc,nchain))
	for(i in 1:nsite){
          for(j in 1:nmcmc){
            for(k in 1:nchain){
              try({
                gsd<-get.stat.dist(P[,,i,j,k])
                w[,i,j,k]<-gsd$w
                dr[i,j,k]<-gsd$dr},silent=T)
            }
          }
        }
	return(list(w=w,dr=dr))
}




ent.p<-function(P,w){
#entropy of a transition probability matrix
#input
#transition probability matrix P
#stationary distribution w
#output
#entropy
	l<-log(P)
	e<-P*l
        e[P==0]<-0 #0 log(0)=0
	h<-colSums(e)
	h<--sum(h*w)
	return(h)
}

get.stat.dist<-function(P){
	# Calculates the stationary distribution of individual transition probability matrices.
	# Used in stat.dist function.
	eig<-eigen(P)
        ind<- which(abs(eig$values-1)==min(abs(eig$values-1))) #find eigenvalue closest to 1
	w<-abs(eig$vectors[,ind]/sum(eig$vectors[,ind])) #pick and scale corresponding eigenvector
        dr<-1/Mod(eig$values[2]) #damping ratio
	return(list(w=w,dr=dr))
}

ent.w<-function(P,w){
# The normalised entropy of a set of transition probability matrices (P).
#input
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transitin probability matrix for site k, MCMC iteration l, chain m
#w is an ncat x nsite x nmcmc x nchain array of stationary distributions, such that w[,k,l,m] is the stationary distribution for site k, MCMC iteration l, chain m
#output  
#H is a nsite X nmcmc x nchain array of entropies; Hr is a nsite X nmcmc x nchain array of normalised entropies. Hprob is an nsite x nmcmc x nchain logical array, with TRUE if H was NaN but there were no NaN values in the corresponding transition probability matrix. Note we get FALSE if H was NaN simply because there were NaN transition probabilities, because this indicates a problem upstream that we can't solve here, or if there were no problems with H.
	ncat<-dim(P)[1]
	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
	H<-array(dim=c(nsite,nmcmc,nchain))
	Hr<-array(dim=c(nsite,nmcmc,nchain))
        Hprob<-array(dim=c(nsite,nmcmc,nchain))
 	Hmax<--log(1/ncat)
	for(i in 1:nsite){
		for(j in 1:nmcmc){
                  for(k in 1:nchain){
			try({
                          H[i,j,k]<-ent.p(P[,,i,j,k],w[,i,j,k])
                          Hr[i,j,k]<-H[i,j,k]/Hmax
                          if(is.na(H[i,j,k])){
                            Hprob[i,j,k]<- !any(is.na(P[,,i,j,k]))#only a (potentially soluble) problem if H is NaN, but the original P matrix has no NaN values
                          }
                          else{
                            Hprob[i,j,k]<-FALSE
                          }
                        },silent=T)
                      }
		}
	}
	return(list(H=H,Hr=Hr,Hprob=Hprob))
}

e.persist.p<-function(P){
#expected persistence time in each state
#input
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transition probability matrix for site k, MCMC iteration l, chain m
#output
#epP is is a ncat X nsite X nmcmc x nchain array of expected persistence time in each state for each site and MCMC iteration
	ncat<-dim(P)[1]
	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
	epP<-array(dim=c(ncat,nsite,nmcmc,nchain))
	for(i in 1:ncat){
		for(j in 1:nsite){
			for(r in 1:nmcmc){
                          for(k in 1:nchain){
				epP[i,j,r,k]<-1/(1-(P[i,i,j,r,k]))
                              }
                        }
                      }
              }
	return(epP)
}

time.absence<-function(P){
# Calculates the average number of time steps until state 1 is reached from each other state.
#input
#P is the transition probability matrix
#output
#expected time to state 1 from each other state
	P<-t(P)				# Transposition of P matrix because Kemeny and Snell define it the other way round.
	Q<-P[-1,-1]			# Identify Q of canonical form (P) assuming state 1 is the absence state.
	I<-diag(dim(Q)[1])
        if(rcond(I-Q)>1e-12){#don't try inverting the fundamental matrix if it's nearly singular
          N<-solve(I-Q)		# Calculate fundamental matrix once Q is known
          t1<-rowSums(N)	# Calculate expected time until absent state is reached given each starting state.
      }
        else{
          t1<-rep(NA,dim(Q)[1])
        }
        return(t1)
}

time.absence.av<-function(P,w){
# Calculates the average number of time steps until state 1 is reached from each other state.
#input
#wscale is the stationary distribution
#P is the transition probability matrix
#output
#t.abs is the expected number of time steps to state 1 from each other state j
#t.abs.av is the expected number of time steps to state 1, from a state sampled at random from the stationary distribution of non-absent states
	t.abs<-time.absence(P)
	wscale<-w[-1]/sum(w[-1])
	t.abs.av<-sum(t.abs*wscale)
	return(list(t.abs=t.abs,t.abs.av=t.abs.av))
}

all.time.absence<-function(P,w){
#Expected times to absence for each site and MCMC rep
#input
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transitin probability matrix for site k, MCMC iteration l, chain m
#w is an ncat x nsite x nmcmc x nchain array of stationary distributions, such that w[,k,l,m] is the stationary distribution for site k, MCMC iteration l, chain m
#output  
# t.abs.all: an ncat(-1) X nsite X nmcmc x nchain array: expected time steps to state 1 from each other state, for each site and mcmc iteration and chain
# t.abs.av.all returns a nsite X nmcmc x nchain array of average time steps until state 1 is reached from some other state, sampled at random from the stationary distribution of non-absent states
	ncat<-dim(P)[1]
	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
	t.abs.all<-array(dim=c(ncat-1,nsite,nmcmc,nchain))
	t.abs.av.all<-array(dim=c(nsite,nmcmc,nchain))
	for(i in 1:nsite){
		for(j in 1:nmcmc){
                  for(k in 1:nchain){
                    try({
                      T.abs<-time.absence.av(P[,,i,j,k],w[,i,j,k])
                        t.abs.all[,i,j,k]<-T.abs$t.abs
                        t.abs.av.all[i,j,k]<-T.abs$t.abs.av},silent=T)
                  }
                }
              }
	return(list(t.abs.all=t.abs.all,t.abs.av.all=t.abs.av.all))
}

is.ergodic<-function(P){
#is a non-negative square matrix ergodic?
  #based on is.ergodic in Appendix S2 of Stott et al 11, MEE 1:242-252
  eig<-eigen(t(P))  
  ind<- which(abs(eig$values-1)==min(abs(eig$values-1))) #find eigenvalue closest to 1
  v<-as.matrix(eig$vectors[,ind])
  av<-as.matrix(abs(v))#dominant left eigenvector (aka reproductive value)
  return(min(av)>0) #ergodic if no zeros in dominant left eigenvector
}

checkergo<-function(P){
#check ergodicity of all P matrices, and of the matrix that results when we make state 1 absorbing
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transitin probability matrix for site k, MCMC iteration l, chain m
#Note: ergoFL is much faster than this, and tells us what we need to know
#output  
# Pergo: an nsite X nmcmc x nchain array: TRUE if the P matrix for this site, chain, and MCMC rep is ergodic, FALSE otherwise
#Qergo: an nsite x nmcmc x chain array: TRUE if the matrix formed by setting rates of leaving state 1 to zero is ergodic, FALSE otherwise

	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
	Pergo<-array(dim=c(nsite,nmcmc,nchain))
	Qergo<-array(dim=c(nsite,nmcmc,nchain))        
	for(i in 1:nsite){
		for(j in 1:nmcmc){
                  for(k in 1:nchain){
                    Pijk<-P[,,i,j,k]
                    try(Pergo[i,j,k]<-is.ergodic(Pijk),silent=T)
                    Pijk[,1]<-0
                    Pijk[1,1]<-1 #state 1 absorbing
                    try(Qergo[i,j,k]<-is.ergodic(Pijk),silent=T)                    
                  }
                }
              }
	return(list(Pergo=Pergo,Qergo=Qergo))
}

ergoFL<-function(P){
#check ergodicity of all P matrices
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transition probability matrix for site k, MCMC iteration l, chain m
#output  
#nsite X nmcmc x nchain arrays: each element is a list of groups of states (FLSL), status of these (FLS), and number of ergodic sets in the full matrix (FLnergo) and the matrix with state 1 set to be absorbing (FLnergoQ)

	nsite<-dim(P)[3]
	nmcmc<-dim(P)[4]
        nchain<-dim(P)[5]
        FLSL<-array(list(NULL),dim=c(nsite,nmcmc,nchain))
        FLS<-array(list(NULL),dim=c(nsite,nmcmc,nchain))
        FLnergo<-array(dim=c(nsite,nmcmc,nchain))
        FLnergoQ <- FLnergo
	for(i in 1:nsite){
		for(j in 1:nmcmc){
                  for(k in 1:nchain){
                    Pijk<-P[,,i,j,k]
                    try({
                      FL<-FoxLandi(Pijk)
                      FLSL[[i,j,k]]<-FL$SL
                      FLS[[i,j,k]]<-FL$S
                      FLnergo[i,j,k]<-FL$nergo
                      #now make state 1 absorbing and repeat the calculation
                      Pijk[,1]<-0
                      Pijk[1,1]<-1 #state 1 absorbing
                      FL<-FoxLandi(Pijk)
                      FLnergoQ[i,j,k]<-FL$nergo
                    },silent=T)
                  }
                }
              }
	return(list(FLSL=FLSL,FLS=FLS,FLnergo=FLnergo,FLnergoQ=FLnergoQ))
}

ergoFL2<-function(P){
#check ergodicity of all P matrices
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transition probability matrix for site k, MCMC iteration l, chain m
#output  
#nsite X nmcmc x nchain arrays: each element is a list of groups of states (FLSL), status of these (FLS), and number of ergodic sets in the full matrix (FLnergo) and the matrix with state 1 set to be absorbing (FLnergoQ)
#this version: builds a list of connection patterns to avoid doing Fox-Landi too many times. Returned as clistindex. Patterns are coded as integers.

  nstate<-dim(P)[1]
  nsite<-dim(P)[3]
  nmcmc<-dim(P)[4]
  nchain<-dim(P)[5]
  bbin<-2^(0:(nstate^2-1))
  FLSL<-array(list(NULL),dim=c(nsite,nmcmc,nchain))
  FLS<-array(list(NULL),dim=c(nsite,nmcmc,nchain))
  FLnergo<-array(dim=c(nsite,nmcmc,nchain))
  FLnergoQ <- FLnergo
  clistindex<-list() #will fill this with unique connection patterns
  clistFL<-list()
  for(i in 1:nsite){
    for(j in 1:nmcmc){
      for(k in 1:nchain){
         Pijk<-P[,,i,j,k]
        if(!any(is.na(Pijk))){
          indB<-codepattern(Pijk,bbin)
          listind<-which(clistindex==indB)
          if(length(listind)>0){#this connection pattern already seen
            FLSL[[i,j,k]]<-clistFL[[listind]]$FLSL
            FLS[[i,j,k]]<-clistFL[[listind]]$FLS
            FLnergo[i,j,k]<-clistFL[[listind]]$FLnergo
            FLnergoQ[i,j,k]<-clistFL[[listind]]$FLnergoQ          
          }else
          {
#            try({
              FL<-NULL
              FL2<-NULL
              FL<-FoxLandi(Pijk)
              FLSL[[i,j,k]]<-FL$SL
              FLS[[i,j,k]]<-FL$S
              FLnergo[i,j,k]<-FL$nergo
              Pijk[,1]<-0  #now make state 1 absorbing and repeat the calculation
              Pijk[1,1]<-1 #state 1 absorbing
              FL2<-FoxLandi(Pijk)
              FLnergoQ[i,j,k]<-FL2$nergo
#            },silent=T)
            ll<-length(clistindex)
            clistindex[ll+1]<-indB
            clistFL[[ll+1]]<-list(FLSL=FL$SL,FLS=FL$S,FLnergo=FL$nergo,FLnergoQ=FL2$nergo)
          }
        }
      }
    }
  }
  return(list(FLSL=FLSL,FLS=FLS,FLnergo=FLnergo,FLnergoQ=FLnergoQ,clistindex=clistindex))
}


FoxLandi<-function(P){
  #Fox-Landi algorithm to identify ergodic sets in a transition probability matrix
  #input P, an n x n transition probability matrix (columns are source states, rows are destination states: note this is transpose of definition used in paper). Will NO LONGER work if we input B, with elements 1 if p_{ij}>0, 0 otherwise (because it now checks that we have a stochastic matrix).
  #outputs a list containing S (each element is a list of states that form a group), SL (each group labelled as Ergodic or Transient), and nergo (number of ergodic sets)
  #source: Fox and Landi (1968) Communications of the ACM 11:619-621

  if(any(is.na(P)) | any(colSums(P)<1-1e-9) | any(P<0)){#we don't know the communication matrix
    return(list(S=list(NULL),SL=list(NULL),nergo=NaN))
  }
  B<-P>0
  n<-dim(P)[1]
  if(!any(B==0)){#first, does every state communicate with every other? If so, quick exit
    SL<-as.list('E') #only one ergodic set
    S<-list(1:n) #set contains all states
    FL<-list(S=S,SL=SL,nergo=getnergo(SL))
  }
  else{
    FL<-FLmain(B,n)
  }
  return(FL)
}

getnergo<-function(SL){#how many ergodic sets in list?
  return(sum(SL=="E"))
}

FLmain<-function(B,n){#does most of the Fox-Landi work
  S<-as.list(1:n)
  SL<-as.list(rep(NA,n))
  abss<-rulea(B)   #find absorbing states
  SL[abss$wa]<-'E' #absorbing states are ergodic
  transs<-ruleb(B,abss$wa,n) #find transient states that communicate with these
  SL[transs$wt]<-'T'
  r<-0 #end of step 1
  while(any(is.na(SL))){#step 2, else if not all states labelled
    if(r==0){#step 3
      i<-which(is.na(SL))[1] #begin a chain with first unlabelled state
      chain<-i
    }
    Bcol<-B[,i]
    Bcol[i]<-FALSE #don't count self link
    irplus<-which(Bcol)[1] #state with which i communicates
    chain<-c(chain,irplus) #extend chain
    if(SL[irplus]=="T"){
      SL[chain]<-"T" #by rule c
      r<-0
    }
    else if(is.na(SL[irplus]) & sum(chain[1:(length(chain)-1)]==irplus)>0){#step 4
      k<-which(chain[1:(length(chain)-1)]==irplus)
      B[,irplus]<-rowSums(B[,chain])>0 #logical OR on cols in chain
      B[irplus,]<-colSums(B[chain,])>0 #logical OR on rows in chain
      B<-B[-chain[(k+1):(length(chain)-1)],-chain[(k+1):(length(chain)-1)]] #delete rows and cols corresponding to elements k+1,k+2,...irplus of chain
      SL<-SL[-chain[(k+1):(length(chain)-1)]]
      S[[chain[k]]]<-unlist(S[chain[k:(length(chain)-1)]])#union of S for items k,k+1,...
      S<-S[-chain[(k+1):(length(chain)-1)]]#this loses all but first element

      absk<-rulea(B)
      if(absk$absorb[irplus]){#state k absorbing
        if(is.matrix(B)){#if only one state left, no need to do this
          trk<-ruleb(B,irplus,dim(B)[1])
          SL[trk$wt]<-"T" #states communicating with this set transient
        }
        SL[irplus]<-"E"
        r<-0
      }
      else{
        r<-k
        chain<-chain[1:k]
        i<-irplus
      }
    }
    else{
      r<-r+1
      i<-irplus
    }
  }
  return(list(S=S,SL=SL,nergo=getnergo(SL)))
}

rulea<-function(B){
  #identify absorbing states
  #input a Boolean matrix with 1 if there is a nonzero transition probability from state j to state i, and 0 otherwise
  #output a Boolean vector absorb with TRUE if a state is absorbing and FALSE otherwise, and the indices wa of absorbing states
  if(!is.matrix(B)){#only one element left
    B<-matrix(B,nrow=1,ncol=1)
  }
  absorb<-(colSums(B)==1 & diag(B)==1)
  wa<-which(absorb)
  return(list(absorb=absorb,wa=wa))
}

ruleb<-function(B,absind,n){
  #identify transient states given index of an absorbing state
    #input a Boolean matrix B with 1 if there is a nonzero transition probability from state j to state i, and 0 otherwise, the indices absind of absorbing states, and the number of states n
  #output a Boolean vector transient with TRUE if a state is transient and FALSE otherwise, and the indices wt of transient states
  Bns<-B
  diag(Bns)<-FALSE #exclude self transitions
  Bns<-Bns[absind,] #absorbing states
  tsum<-colSums(matrix(Bns,ncol=n)) #communicate directly with any absorbing state
  transient<-tsum >0
  wt<-which(transient)
  return(list(transient=transient,wt=wt))
}

testB<-function(){#make some test cases
  #communication matrix in the Fox/Landi paper)
  B<-matrix(c(c(0,0,0,1,0,0,0,0,0),c(0,1,0,0,0,1,1,0,0),c(0,1,0,0,0,0,0,1,0),c(1,0,0,0,0,0,0,0,1),c(0,1,0,0,1,0,0,0,0),c(0,0,0,0,0,0,1,0,0),c(0,0,0,0,0,0,1,0,0),c(0,0,1,0,0,0,0,0,0),c(1,1,0,0,0,0,0,0,1)),nrow=9)
  B1<-t(B) #transpose for population biology
  #correct output: {1,4,9},{3,8},{5} ergodic, 2, 6, and 7 transient

  #communication matrix where all states communicate
  B2<-matrix(nrow=9,ncol=9)
  B2[,]<-1
  #correct output: {1,2,...,9} ergodic

  #one pair of states has no direct communication
  B3<-matrix(1,nrow=4,ncol=4)
  B3[4,1]<-0
  #correct output: {1,2,...,4} ergodic

  P1<-matrix(c(0,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0),nrow=4) #correct output: {1,2},{3,4} ergodic
  P2<-matrix(c(0.8,0.2,0,0,0,0.8,0.2,0,0.2,0,00.8,0,0,0,0,1),nrow=4) #correct output: {1,2,3},{4} ergodic
  P3<-matrix(c(0.8,0,0.2,0,1,0,0,0,0.2,0,0.8,0,0,0,0,1),nrow=4) #correct output: {1,3},{4} ergodic, 2 transient
  P4<-matrix(c(0.8,0,0.2,0,0.9,0,0,0.1,0.2,0,0.8,0,0,0.1,0,0.9),nrow=4) #correct output: {1,3} ergodic, 2 and 4 transient
  return(list(B1=B1,B2=B2,B3=B3,P1=P1,P2=P2,P3=P3,P4=P4))
}

checkdef<-function(P){
#check whether there are any undefined elements in all P matrices (and that they have column sums (close to) 1, and no negative elements
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transition probability matrix for site k, MCMC iteration l, chain m
#output a list:
#Pdef is nsite x nmcmc x nchain: TRUE if all elements are defined, FALSE otherwise
#Pvalid is nsite x nmcmc x nchain: TRUE if we have a stochastic matrix (col sums 1, all elements non-negative), FALSE otherwise (or NA if Pdef is false).

	nsite <- dim(P)[3]
	nmcmc <- dim(P)[4]
        nchain <- dim(P)[5]
	Pdef <- array(FALSE,dim=c(nsite,nmcmc,nchain))
        Pvalid <- Pdef
	for(i in 1:nsite){
		for(j in 1:nmcmc){
                  for(k in 1:nchain){
                    Pijk <- P[,,i,j,k]
#                    try(Pdef[[i,j,k]] <- !any(is.na(Pijk)),silent=T)
                    Pdef[[i,j,k]] <- !any(is.na(Pijk))                    
                    Pvalid[[i,j,k]] <- !any(colSums(Pijk)<1-1e-9) & !any(Pijk<0)
                  }
                }
              }
	return(list(Pdef=Pdef,Pvalid=Pvalid))
}

countpatterns<-function(P){
#how many different connection patterns?
#P is an ncat x ncat x nsite x nmcmc x nchain array, such that P[,,k,l,m] is the transition probability matrix for site k, MCMC iteration l, chain m
#output a list containing uindB, the unique connection patterns (represented as integers), and indB, an nsite x nmcmc x nchain array of connection patterns
  nstate<-dim(P)[1]
  nsite<-dim(P)[3]
  nmcmc<-dim(P)[4]
  nchain<-dim(P)[5]
  indB<-array(dim=c(nsite,nmcmc,nchain))
  bbin<-2^(0:(nstate^2-1))
  for(i in 1:nsite){
    for(j in 1:nmcmc){
      for(k in 1:nchain){
        Pijk<-P[,,i,j,k]
        B<-as.integer(Pijk>0)
        indB[i,j,k]<-sum(B*bbin) #turn the connection pattern into an integer code
      }
    }
  }
  uindB<-as.numeric(names(table(indB))) #extract the unique connection patterns
  return(list(uindB=uindB,indB=indB))
}

#input a transition probability matrix P and the vector
#  bbin<-2^(0:(nstate^2-1))
#return an integer code identifying the pattern of nonzero transition probabilities
codepattern<-function(P,bbin){
  B<-as.integer(P>0)
  indB<-sum(B*bbin) #turn the connection pattern into an integer code
  return(indB)
}
