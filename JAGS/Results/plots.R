rm(list=ls())
setwd("C:/Users/Robert/Documents/PhD/R files/DS_model_weeds/Bayesian model/Results")# set path to the distribution folder

trans.rotation_09<-read.table(file="rot09.txt",sep="\t",dec=".")
trans.norotation_09<-read.table(file="norot09.txt",sep="\t",dec=".")
naive<-read.table(file="naive09.txt",sep="\t",dec=".")

trans.rotation_08<-read.table(file="rot08.txt",sep="\t",dec=".")
trans.norotation_08<-read.table(file="norot08.txt",sep="\t",dec=".")


trans.rotation_09<-as.matrix(trans.rotation1)
trans.norotation_09<-as.matrix(trans.norotation1)
trans.rotation_08<-as.matrix(trans.rotation_08)
trans.norotation_08<-as.matrix(trans.norotation_08)
trans.naive<-as.matrix(naive)


par(mfrow=c(3,2))
par(mar=c(5, 4, 4, 2))

image(z=t(trans.rotation_09), axes=F, col=rainbow(n=40,start=0.5, end=0.6),xlab="Density state in 2008",ylab="Density state in 2009",main="Rotation")
mtext(side=1,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)
mtext(side=2,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)

image(z=t(trans.norotation_09), axes=F, col=rainbow(n=40,start=0.5, end=0.6),xlab="Density state in 2008",ylab="Density state in 2009",main="No Rotation")
mtext(side=1,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)
mtext(side=2,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)

image(z=t(trans.rotation_08), axes=F, col=rainbow(n=40,start=0.5, end=0.6),xlab="Density state in 2008",ylab="Density state in 2009",main="Rotation")
mtext(side=1,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)
mtext(side=2,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)

image(z=t(trans.norotation_08), axes=F, col=rainbow(n=40,start=0.5, end=0.6),xlab="Density state in 2008",ylab="Density state in 2009",main="No Rotation")
mtext(side=1,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)
mtext(side=2,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)

image(z=t(trans.naive), axes=F, col=rainbow(n=40,start=0.5, end=0.6),xlab="Density state in 2008",ylab="Density state in 2009",main="No Rotation")
mtext(side=1,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)
mtext(side=2,c("Abs","L","M","H","V"),at=c(0,0.25,0.5,0.75,1),line=0.6)


pred1<-trans.rotation
pred2<-trans.norotation
col.vect<-gray(0:8/8)

par(mfrow=c(1,1))
plot(sqrt(trans.rotation[1,]),sqrt(trans.norotation[1,]),pch=16,xlab="sqrt transitions when no rotation",ylab="sqrt transitions when rotation",col=col.vect[1],xlim=c(0,1),ylim=c(0,1),mgp=c(1.8,0.6,0),las=1,main="Predicted rotation effect 08-09",bty="n", yaxs='i', xaxs='i')
points(sqrt(trans.rotation[2,]),sqrt(trans.norotation[2,]),pch=16,col=col.vect[5])
points(sqrt(trans.rotation[3,]),sqrt(trans.norotation[3,]),pch=16,col=col.vect[6])
points(sqrt(trans.rotation[4,]),sqrt(trans.norotation[4,]),pch=16,col=col.vect[7])
points(sqrt(trans.rotation[5,]),sqrt(trans.norotation[5,]),pch=16,col=col.vect[8])
points(sqrt(trans.rotation[2,]),sqrt(trans.norotation[2,]),pch=1)
points(sqrt(trans.rotation[3,]),sqrt(trans.norotation[3,]),pch=1)
points(sqrt(trans.rotation[4,]),sqrt(trans.norotation[4,]),pch=1)
points(sqrt(trans.rotation[5,]),sqrt(trans.norotation[5,]),pch=1)
abline(0,1)


legend("topleft",col=col.vect[-c(2,3,4)],pch=16,legend=(c("to Abs","to L","to M","to H","to V")),bty="n")
legend("topleft",col="black",pch=1,legend=(c("","","","","")),bty="n")

# function to work out equilibrium distributions for average matrices. 
state.distrib<-function(av.mat){ 

  av.mat<-t(av.mat)
    # initials states
  v<-c(1,0,0,0,0)
  
  # second vector with arbitrary values which get overwritten in loop below
  vec2<-c(1,1,1,1,1)
  
  #first iteration
  vec1<-v%*%av.mat
  
  
  #counter and 5 density states
  counter=1
  states<-cbind("1","2","3","4","5")
  
  # loop that runs iterations until equilibrium is reached - ie no change between vh density state between iterations
  while ((vec1[5]/vec2[5])!=1){
    counter<-counter+1;
    vec2 <- vec1%*%av.mat;
    vec1 <- vec2%*%av.mat;
    if (counter==200) break
    
  }
  
  # combine vectors into data frame and plot
  f<-as.vector(vec1)
  s<-as.vector(states)
  final<-data.frame(s,f)
  barplot(final$f,names.arg=c("A","L","M","H","VH"), ylim=c(0,1), ylab=("Proportion"),xlab="Density state")
}



state.distrib(trans.rotation_09)
state.distrib(trans.norotation_09)
state.distrib(trans.rotation_08)
state.distrib(trans.norotation_08)
state.distrib(naive)


damp.ratio<- function(av.matrix){
  eigen.vals<- eigen(av.matrix)
  DR<-Mod(eigen.vals$values)
  DR2 <- DR[1]/abs(DR[2])
  
}

damping_ratios <- lapply(())

