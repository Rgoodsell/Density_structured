transition.count<-function(spdata,ncat){
#count the number of transitions between each state
#input
#spdata is a data frame of MarClim data
#ncat is number of states
#output
#m is a matrix of transition counts (cols are source, rows are destination)
m<-matrix(0,nrow=ncat,ncol=ncat)
	for(n in 1:dim(spdata)[1]){
		i<-spdata$ytplus[n]
		j<-spdata$yt[n]
		m[i,j]<-m[i,j]+1
	}
return(m)
}
