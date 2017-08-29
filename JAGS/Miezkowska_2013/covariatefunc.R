
fmfa <- multinom( year2008 ~ year2007 + farm, data=dat_rotation)
fmr <- multinom( year2008 ~ year2007 + rotation,data=dat_rotation)


fmfar<- multinom( year2008 ~ year2007 + rotation + farm, data=dat_rotation)
fmfi <- multinom (year2008 ~ year2007 + field, data=dat_rotation, MaxNWts= 2361)
fmfir<- multinom(year2008 ~ year2007 + field +rotation, data=dat_rotation, MaxNWts = 2366)

fm.out <- trans.mat.cov(fm, dat_rotation$farm)

trans.mat <- function(vobj) {
  # Extract the predictors
  lp1 <- c(0,coef(fm)[1:4])
  lp2 <- coef(fm)[5:20]
  
  # Set up as a matrix of interaction coefficients
  lp2 <- matrix(lp2, nrow =4 , ncol = 4, byrow = FALSE)
  lp2 <- cbind( rep(0, 4), lp2)
  lp2 <- rbind( rep(0, 5),lp2)
  #that creates a matrix with the interaction plus an added 5th column and 5th row with zeros
  
  p <- matrix(0, nrow = 5, ncol = 5)
  
  # Make the predictions and apply
  for(i in 1:5){ 
    pr <- lp1+  lp2[,i] 
    p[i,] <- exp(pr) / sum(exp( pr ))
  }	
  return(p)
}



trans.mat.cov <- function(vobj, x) {
  
  list.res<-vector("list",length(levels(x)))
  
  # Extract the predictors
  # vglm takes the last level of the predictor as the reference
  lp1 <- c(0,coef(vobj)[1:4])
  # this differs from vglm
  lp2 <- coef(vobj)[5:20]
  lp3 <- coef(vobj)[21:length(coef(vobj))]
  
  # Set up as a matrix of interaction coefficients
  lp2 <- matrix(lp2, nrow =4 , ncol = 4, byrow = FALSE)
  lp2 <- cbind( rep(0, 4), lp2)
  lp2 <- rbind(rep(0, 5),lp2)
  # this last line differs from vglm
  
  # The covariate slopes
  lp3 <-matrix(lp3, nrow =4 , ncol = length(lp3)/4, byrow = FALSE)
  lp3 <- cbind( rep(0, 4), lp3)
  lp3 <- rbind( rep(0,ncol(lp3)),lp3)
  
  # Make the predictions and apply
  counter<-1
  for(j in 1:length(levels(x))){
    p <- matrix(0, nrow = 5, ncol = 5)
    for(i in 1:5){ 
      pr <- lp1+  lp2[,i] + lp3 [,j]
      p[i,] <- exp(pr) / sum(exp( pr ))
    }	
    list.res[[counter]]<-p
    counter<-counter+1
  }
  names(list.res)<-levels(x)
  return(list.res)
}
