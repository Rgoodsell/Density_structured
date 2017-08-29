#filled contour plot of three variables
#input x and y, the x and y coordinates as vectors (same length)
#z is a vector of values of the z variable at pairs x,y
#xlab and ylab are axis labels
envelope.plot <- function(x,y,z,xlab,ylab)
{
	nstep <- length(x)
	filled.contour(x=x,y=y,z=matrix(z,nrow=nstep,ncol=nstep,byrow=FALSE),xlab=xlab,ylab=ylab,col=gray.colors(19),nlevels=20) 
}

#plot mean and standard deviation of a function of two variables
#input x and y, the x and y coordinates as vectors (same length)
#zmean is a vector of means at the specified coordinates
#zsd is a vector of standard deviations at the specified coordinates
#xlab and ylab are x and y axis labels
#contour: set to TRUE to draw contour lines
#region: set to TRUE to shade to represent z value
#maina and mainb are panel labels for the two panels
#if fixat is TRUE (default), we draw contours at 0.05-intervals between 0 and 1 (otherwise let contourplot decide where to put them
envelope.plot2 <- function(x,y,zmean,zsd,xlab,ylab,contour=FALSE,region=FALSE,maina,mainb,fixat=TRUE)
{
	if(fixat==TRUE){
		plot1<-contourplot(zmean ~ x*y,contour=contour,region=region,at=seq(from=0,to=1,by=0.05),col.regions=grey.colors,xlab=xlab,ylab=ylab,main=list(maina,cex=1,font=1),scales=list(tck=-0.5))
		plot2<-contourplot(zsd ~ x*y,contour=contour,region=region,at=seq(from=0,to=1,by=0.05),col.regions=grey.colors,xlab=xlab,ylab=ylab,main=list(mainb,cex=1,font=1),scales=list(tck=-0.5))
	} else {
		plot1<-contourplot(zmean ~ x*y,contour=contour,region=region,col.regions=grey.colors,xlab=xlab,ylab=ylab,main=list(maina,cex=1,font=1),scales=list(tck=-0.5))
		plot2<-contourplot(zsd ~ x*y,contour=contour,region=region,col.regions=grey.colors,xlab=xlab,ylab=ylab,main=list(mainb,cex=1,font=1),scales=list(tck=-0.5))

}
	print(plot1,position=c(0, .5, 1, 1), more=TRUE)
	print(plot2,position=c(0, 0, 1, .5))
}

#ENDED UP NOT USING THIS FUNCTION: STILL NEEDS A LOT OF WORK
#contour panel function
panel.contour <- function(x,y,z,contour,region,cxlabi,subscripts=1,...){
#less space around each panel
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))

         	panel.contourplot(x,y,z,contour=contour,region=region,subscripts=subscripts,at=seq(from=0,to=1,by=0.1),col.regions=grey.colors,scales = list(tck=-0.5,x=list(cex=0.75),y=list(cex=0.75)),par.settings=theme.novpadding,xlab=list(xlab,cex=0.75,col=cxlabi),ylab=list(ylab,cex=0.75))


}

#plot mean and standard deviation of a set of four stationary probabilities
#and their standard deviations
#as functions of two variables
#input x and y, the x and y coordinates as vectors (same length)
#zmeans is list of vectors means at the specified coordinates for each stationary prob
#zsds is list of vectors of standard deviations at the specified coordinates for each stationary prob
#xlab and ylab are x and y axis labels
#contour: set to TRUE to draw contour lines
#region: set to TRUE to shade to represent z value
#mains is list of panel labels (order mean than standard deviation, states 1 to 4)
statcontour.plot <- function(x,y,zmeans,zsds,xlab,ylab,contour=FALSE,region=FALSE,mains)
{
	plots <- vector("list",8)

#less space around each panel
theme.novpadding <-
   list(layout.heights =
        list(top.padding = 0,
 	    main.key.padding = 0,
 	    key.axis.padding = 0,
 	    axis.xlab.padding = 0,
 	    xlab.key.padding = 0,
 	    key.sub.padding = 0,
 	    bottom.padding = 0),
        layout.widths =
        list(left.padding = 0,
 	    key.ylab.padding = 0,
 	    ylab.axis.padding = 0,
 	    axis.key.padding = 0,
 	    right.padding = 0))
	for(i in 1:4){#trying to force all panels to have same height of axes
	      if(i==4){
		cxlabi="black"} else {
		cxlabi="transparent"
		}

#CALL TO PANEL.CONTOUR NOT WORKING YET
#		plots[[2*i-1]] <- panel.contour(x=z ~ x*y, data = data.frame(x=as.matrix(x),y=as.matrix(y),z=as.matrix(zmeans[[i]])),contour=contour,region=region,cxlabi)

		plots[[2*i-1]]<-contourplot(zmeans[[i]] ~ x*y,contour=contour,region=region,at=seq(from=0,to=1,by=0.1),col.regions=grey.colors,scales = list(tck=-0.5,x=list(cex=0.75),y=list(cex=0.75)),par.settings=theme.novpadding,xlab=list(xlab,cex=0.75,col=cxlabi),ylab=list(ylab,cex=0.75),main=list(label=mains[[2*i-1]],cex=0.75))

		plots[[2*i]]<-contourplot(zsds[[i]] ~ x*y,contour=contour,region=region,at=seq(from=0,to=1,by=0.1),col.regions=grey.colors,scales = list(tck=-0.5,x=list(cex=0.75),y=list(cex=0.75)),par.settings=theme.novpadding,xlab=list(xlab,cex=0.75,col=cxlabi),ylab=list(ylab,cex=0.75,col="transparent"),main=list(label=mains[[2*i]],cex=0.75))
	}

	print(plots[[1]],position=c(0, .75, .5, 1), more=TRUE)
	print(plots[[2]],position=c(.5, .75, 1, 1), more=TRUE)
	print(plots[[3]],position=c(0, .5, .5, .75), more=TRUE)
	print(plots[[4]],position=c(.5, .5, 1, .75), more=TRUE)
	print(plots[[5]],position=c(0, .25, .5, .5), more=TRUE)
	print(plots[[6]],position=c(.5, .25, 1, .5), more=TRUE)
	print(plots[[7]],position=c(0, 0, .5, .25), more=TRUE)
	print(plots[[8]],position=c(.5, 0, 1, .25))
}



#calculate values for plotting relationship between a univariate summary statistic and two explanatory variables
#inputs
#w is sites x reps x chains array of summary statistic
#cd is sites x reps x chains array of TRUE/FALSE values: are all elements in corresponding transition probability matrix defined?
#FL is sites x reps x chains array of output from ergoFL, includes FL$nergo, number of ergodic sets in transition probability matrix
#output a list containing the following vectors, each with one element per site:
#smean, mean summary statistic (over reps and chains)
#ssd, standard deviation of summary statistic (over reps and chains)
#pdiscard: proportion of reps we threw away (not all values defined, or more than one ergodic set)

envelope.plotvalues<-function(w,cd,FL){
  nsite <- dim(w)[1]
  smean <- array(dim=c(1,nsite))
  ssd <- array(dim=c(1,nsite))
  pdiscard<-0
  for(i in 1:nsite){
    wi <- w[i,,]
    cdi <- cd[i,,]
    FLi <- FL$FLnergo[i,,]==1
    wsub <- wi[cdi & FLi]#subset: all elements defined and one ergodic set
    smean[i] <- mean(wsub)
    ssd[i] <- sd(wsub)
    pdiscard <- pdiscard + length(wi) - length(wsub)  #record proportion of reps we threw away!
  }
  pdiscard<-pdiscard/prod(dim(w))
  return(list(smean=smean,ssd=ssd,pdiscard=pdiscard))
}

