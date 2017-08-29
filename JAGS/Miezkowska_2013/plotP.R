plotP<-function(ncat,P,iters,spdata,sourcecat){
  #marginal plots of predicted transition probabilities against explanatory variables
  #inputs
  #ncat is number of abundance categories
  #P is an array of transition probability matrices, with dimensions destination state, source state, site, MCMC iteration,chain. Thus P[,,2,7,3] is the transition probability matrix for the second site, 7th MCMC iteration, 3rd chain
  #iters is the iterations we want to average over (e.g.1:100)
  #spdata is the raw data
  #sourcecat is the source category
  #outputs
  #scatter plots of mean transition probabilities against explanatory variables. Each circle is the posterior mean transition probability (over chains and iterations) from sourcecat to each other category i for a given row of the raw data (a site in a particular pair of years). + signs are at 1 if in a given row, there was a transition from sourcecat to i, and 0 otherwise. The line is the prediction from a binomial GAM (used as a smoother, fitted separately for each destination category)
  #note, fetch is divided by 5 throughout (we realized that it's not really in km, it's in 200m grid cells)

  par(mfrow=c(ncat,2))
  par(mar=c(5,4,2,2))
  subdata<-spdata[spdata$yt==sourcecat,] #subset of raw data in which the source state is the one of interest
  Psub<-apply(P[,,spdata$yt==sourcecat,iters,],c(1,2,3),mean,na.rm=TRUE)[,sourcecat,] #mean transition probabilities out of the source state of interest, for each row of the 
#  panellab<-c("a","b","c","d","e","f","g","h")
  panellab<-letters[1:I(2*ncat)]
  for(i in 1:ncat){

    mod <- gam( (ytplus==i)  ~ s( sst ), data = subdata, family = binomial ) #binomial gam (needs mgcv library)
    plabel<-paste(as.character(i),as.character(sourcecat),sep="")
    par(tcl=0.1)#tick marks inwards and short
    plot(subdata$sst,Psub[i,],ylim=c(0,1),xlab=expression(paste("Mean winter SST (",degree,"C)")),ylab=substitute(italic(p)[plabel],list(plabel=plabel)))
 #   points(subdata$sst,subdata$ytplus==i,pch=3) #but points get hidden by overplotting
    phist <- plotPhist(subdata$sst,subdata$ytplus,i)    #make little histograms at each value of explanatory variable
    points(subdata$sst,phist$yy0,pch=20,cex=0.2)
    points(subdata$sst,phist$yy1,pch=20,cex=0.2)

    newdata<-data.frame(sst=seq(from=min(subdata$sst),to=max(subdata$sst),length.out=101))
    gampred<-predict(mod,newdata,type="response")
    lines(newdata$sst,gampred)
    axl<-par("usr") #get axis limits
    par(xpd=NA) #allow text outside limits
    text(axl[1]-0.1*(axl[2]-axl[1]),axl[4]+0.1*(axl[4]-axl[3]),panellab[2*i-1])

    mod <- gam( (ytplus==i)  ~ s( fetch), data = subdata, family = binomial ) #binomial gam (needs mgcv library)
    plabel<-paste(as.character(i),as.character(sourcecat),sep="")
    par(tcl=0.1)
    plot(subdata$fetch/5,Psub[i,],ylim=c(0,1),xlab="Wave fetch (km)",ylab=substitute(italic(p)[plabel],list(plabel=plabel)))
#    points(subdata$fetch/5,subdata$ytplus==i,pch=3)
    phist <- plotPhist(subdata$fetch/5,subdata$ytplus,i)    #make little histograms at each value of explanatory variable
    points(subdata$fetch/5,phist$yy0,pch=20,cex=0.1)
    points(subdata$fetch/5,phist$yy1,pch=20,cex=0.1)

    newdata<-data.frame(fetch=seq(from=min(subdata$fetch),to=max(subdata$fetch),length.out=101))
    gampred<-predict(mod,newdata,type="response")
    lines(newdata$fetch/5,gampred)
        axl<-par("usr") #get axis limits
    par(xpd=NA) #allow text outside limits
    text(axl[1]-0.1*(axl[2]-axl[1]),axl[4]+0.1*(axl[4]-axl[3]),panellab[2*i])

  }
}

#data for plotting little histograms of response subresp (1 if we see a specified transition (given by i, the destination state), 0 otherwise) at each value of an explanatory variable subexp.
plotPhist<-function(subexp,subresp,i){

    np <- dim(subexp)[1]
    yy0 <- array(NA,dim=c(np,1))
    yy1 <- yy0
    step <- 0.02
    stepstart <- 0.07
    yy0[!(subresp==i)] <- -stepstart
    yy1[subresp==i] <- 1+stepstart
    for(j in unique(subexp)){
      npj0 <- sum(subexp==j & !(subresp==i))
      npj1 <- sum(subexp==j & subresp==i)      

      if(npj0>0) yy0[subexp==j & !(subresp==i)] <- seq(from=-stepstart,to= -stepstart - (npj0-1)*step,by= -step)
      if(npj1>0) yy1[subexp==j & subresp==i] <- seq(from=1+stepstart,to= 1+stepstart+(npj1-1)*step,by= step)
    }
    return(list(yy0=yy0,yy1=yy1))
}
