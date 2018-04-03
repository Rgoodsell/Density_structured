########### Simulate & compare densities from fitted gamma distribution ###########
# Functions to integrate gamma dist,optimise, fit, & plot #


##### Integrates gamma function #####
gamma.int.fun <- function( breaks, shape, rate ) {
  
  # breaks = the break points that define the numerical abundances contained within each density state
  # shape  = the shape parameter of the gamma distribution to be fitted to the stable density structure
  # rate   = the rate parameter of the gamma distribution to be fitted to the stable density structure
  # both rate and shape are determined using optim below
  
  n <- length(breaks) # Number of density states
  starts <- breaks[1:(n-1)] 
  ends <-  breaks[ 2:n]  
  probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate) # calculates definite integral for a set of shape & rate parameters 
  return( probs )
}

##### Optimizes the gamma function for data and given breaks #####
opt.gamma.fun <- function( ds, breaks, init = c(1, 1) ) {
  
  # --- Function that is optimised, crudely based on least squares --- #
  optFun <- function( pars) {
    
    # takes values for shape & rate parameters of gamma distribution
    # returns sums of squares value  - to be minimised using optim below 
    
    shape <- pars[1]
    rate <-pars[2]
    probs <- gamma.int.fun( breaks, shape, rate)
    ssq <-  ( ds - probs )^2 #  residuals squared 
    
    return( sum( ssq) )
  }
  
  ds <- ds / sum(ds) # proportion of counts in each density category 
  op <- optim( init, optFun) # finds the best fit of shape and rate parameters
  n <- length(breaks) 
  shape <- op$par[1] 
  rate <- op$par[2]
  starts <-  breaks[1:(n-1)] 
  ends <-  breaks[ 2:n] 
  
  probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate) # integral  
  
  
  ## Code to plot 
  
  # x <- seq(min(breaks), max(breaks), by = (max(breaks) - min(breaks)) / 100)
  # 
  # fx <-  pgamma(x[2:length(x) ], shape, rate) - pgamma(x[1:(length(x)-1)], shape, rate)
  # 
  # plot( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, probs, type = "l" , log = "x",
  #       xlab = "State", ylab = "Proportion", ylim = c(0,1))
  # 
  # points( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, ds, col= "red" )
  
  mean <- shape / rate # Average density per sample quadrat?
  var <- shape / rate^2 # Variance
  
  
  return(data.frame(mean=mean,var=var))
  #return( list(shape = shape, rate = rate, mean = mean, var = var, probs = probs, ds = ds) )
  
}

#get DS distributions for each field

fitGamma <- function(P,mydata,model,breaks){
  
  # Get a simulated & observed set of DS distributions for each field
  # get number of fields & field level ID 
  nField  <- length(unique(mydata$field))
  raneffID <-  as.numeric(mydata$field)
  K <- length(unique(mydata$ytplus))
  
  # set up output objects
  ds.diffs <- array(0,dim=c(nField,5,2)) # df of differences for each field [f,]
  
  
  
  # Loop through fields
  for(f in 1:nField){
    
    
    
    
    fsub_preds <-  P[raneffID==f,,] # Select theta estimates by random effect index
    fsub_obs <- mydata[raneffID==f,] # get observed distribution of DS in each field
    
    
    obsMat <- data.matrix(table(fsub_obs[,2])) # get observed outcome & source dist
    obs <- obsMat/sum(obsMat) # calculate as proportion
    
    ds.diffs[f,,1] <-  apply(fsub_preds,2,mean)
    ds.diffs[f,1:length(obs),2] <- obs    
    
    
    
  }
  
  
  pred.fits  <- apply(ds.diffs[,,1],1,opt.gamma.fun,breaks=breaks,init=c(1,1)) %>% do.call(rbind,.)
  
  obs.fits  <- apply(ds.diffs[,,2],1,opt.gamma.fun,breaks=breaks,init=c(1,1)) %>% do.call(rbind,.)
  
  
  diffs <- pred.fits[,1] - obs.fits[,1]
  
  return(data.frame(diffs,model=model))
  
}

# Plots the above by model
breaks <- c(0,1,160,450,1450,2000)

plot_gamma_fits <- function(allFits,levels){
  
  #allFits = all gamma distribution fits from above functions
  #levels = levels from model variable in allFits, used for reordering 
  
  # Get summary data
  avgFits<-   allFits %>%  group_by(model) %>% summarise(mean = median(diffs),lower = quantile(diffs,probs = 0.10),
                                                         upper = quantile(diffs,probs=0.80))  
  
  
  # Relevel
  allFits$model <- factor(allFits$model, levels = levels)
  avgFits$model <- factor(avgFits$model , levels = levels)
  
  
  
  # Plot 
  ggplot(avgFits,aes(model,mean))+
    geom_hline(data=avgFits,aes(yintercept = 0),size=1,colour="grey",lty=2)+
    geom_jitter(data=allFits,aes(model,diffs,colour=model),alpha=0.25,size=4,width=0.1)+
    geom_errorbar(data=avgFits,aes(ymax=upper,ymin=lower),size=1,width=0)+
    geom_point(data=avgFits,aes(model,mean),size=5,colour="black") +
    theme_classic()+
    theme(axis.title = element_text(size=30,face="bold"),
          axis.text = element_text(size=28),
          legend.text = element_text(size=28),
          legend.title = element_text(size=30),
          legend.key.height=unit(3,"line"),
          strip.text = element_text(size=28))+
    scale_y_continuous(limits = c(-3, 3))+
    labs(x="Model",y="Error" ,colour="Model")
  
  
  
}
