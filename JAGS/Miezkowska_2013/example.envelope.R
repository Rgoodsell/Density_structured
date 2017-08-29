#"climate envelope" for simulated data

rm(list=ls())
load("example.RData")
library(rjags)
rm(list=lsf.str()) #remove all functions to make sure we don't get old versions by mistake

#construct summary statistics
source("construct.P.R")
source("stat.dist.R")
nstep<-10 #we get an nstep x nstep grid of the two explanatory variables
nrep <- NULL #set to null to use all reps
x1new<-seq(from=min(spdata$scalesst),to=max(spdata$scalesst),length.out=nstep) #nstep equally-spaced values on first variable
x2new<-seq(from=min(spdata$scalefetch),to=max(spdata$scalefetch),length.out=nstep) #nstep equally-spaced values on second variable

x1plot<-seq(from=min(spdata$sst),to=max(spdata$sst),length.out=nstep)
x2plot<-seq(from=min(spdata$fetch/5),to=max(spdata$fetch/5),length.out=nstep) #divide by 5 because fetch is in 200m grid cells, but we want to plot in km

xnew <- as.matrix(expand.grid(sst=x1new,fetch=x2new))
xplot <- as.matrix(expand.grid(sst=x1plot,fetch=x2plot))

print("making Pnew")
Pnew<-construct.p(alljags,xnew,nmcmc=nrep) #get transition probability matrices at these values
print("checkdef")
cdef1<-checkdef(Pnew) #are all elements defined, and do we have a stochastic matrix?
cd1<-cdef1$Pdef & cdef1$Pvalid
print("stat.dist")
w1<-stat.dist(Pnew) #stationary distribution and damping ratio
print("ergoFL2")
FL1<-ergoFL2(Pnew) #how many ergodic sets?
print("all.time.absence")
eta1<-all.time.absence(Pnew,w1$w) #expected time to absence
print("e.persist.p")
epa1<-e.persist.p(Pnew) #expected persistence time in each state
print("ent.w")
entro1<-ent.w(Pnew,w1$w) #entropy

library(lattice)
source("envelope.plot.R")
epv<-envelope.plotvalues(w1$w[1,,,],cd1,FL1) #this gets mean prob of absence at each grid value, over reps and chains

envelope.plot2(x=xplot[,1],y=xplot[,2],zmean=1-epv$smean,zsd=epv$ssd,contour=TRUE,region=FALSE,xlab=expression(paste("Mean winter SST (",degree,"C)")),ylab="Wave fetch (km)",maina="a: probability of presence",mainb="b: s.d. (probability of presence)")


#envelope for entropy
entropv<-envelope.plotvalues(entro1$Hr,cd1,FL1)
envelope.plot2(x=xplot[,1],y=xplot[,2],zmean=entropv$smean,zsd=entropv$ssd,contour=TRUE,region=FALSE,xlab=expression(paste("Mean winter SST (",degree,"C)")),ylab="Wave fetch (km)",maina="a: normalized entropy",mainb="b: s.d. (normalized entropy)")
