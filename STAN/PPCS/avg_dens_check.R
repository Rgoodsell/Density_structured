### Function to check observed vs predicted average densities ### 
avgDens.check <- function(P,mydata,model){
  
  ## Takes:
  # P: site level transition probabilities to source states [,1:K,] at site [n,,] for all iterations [,,i]
  # mydata: data used in model to generate P
  # Character string to ID model
  
  
  ## Returns
  # ds.diffs: df with difference between average field level observed & predicted ds'
  
  # get number of fields & field level ID 
  nField  <- length(unique(mydata$field))
  raneffID <-  as.numeric(mydata$field)
  K <- length(unique(mydata$ytplus))
  
  # set up output objects
  ds.diffs <- array(0,dim=c(nField,1)) # df of differences for each field [f,]
  
  
  
  # Loop through fields
  for(f in 1:nField){
    
    obs.dist <- rep(0,5) # observed DS distribution 
    
    
    fsub_preds <-  P[raneffID==f,,] # Select theta estimates by random effect index
    
    
    fsub_obs <- mydata[raneffID==f,] # get observed distribution of DS in each field
    
    fsub_mean <-  apply(fsub_preds,2,mean)
    
    fsub_sim <- mean(sample(1:5,nrow(fsub_obs),prob=fsub_mean,replace = T)) # simulate from posterior & take mean density state
    #fsub_se <- apply(fsub_preds,2,function(x) sqrt(sd(x))/length(x))
    ds.diffs[f,] <- fsub_sim - mean(fsub_obs[,2]) 
    
    
  }
  ds.diffs  <- as.data.frame.matrix(ds.diffs)
  ds.diffs$model <- model
  return(ds.diffs)
  
  
}


### Function to plot average densities #####
plot_avg_dens <- function(allAvg,levels){
  
  
  ## Takes:
  # allAvg : object containing all average densities from all models
  # levels : levels of models, used to reorder plotting. 
  
  
  
  # Normalise and get summary data
  allAvg$V1 <- allAvg$V1/max(abs(allAvg$V1))
  avgSum <-   allAvg %>%  group_by(model) %>% summarise(mean = mean(V1),lower = quantile(V1,probs = 0.10),
                                                        upper = quantile(V1,probs=0.90))  
  
  
  # Relevel models
  allAvg$model <- factor(allAvg$model, levels = levels)
  avgSum$model <- factor(avgSum$model , levels = levels)
  
  
  # Plot
  P <- ggplot(avgSum,aes(model,mean))+
    geom_hline(data=avgSum,aes(yintercept = 0),size=1,colour="grey",lty=2)+
    geom_jitter(data=allAvg,aes(model,V1,colour=model),alpha=0.25,size=4,width=0.1)+
    geom_errorbar(data=avgSum,aes(ymax=upper,ymin=lower),size=1,width=0)+
    geom_point(size=5,colour="black") +
    theme_classic()+
    theme(axis.title = element_text(size=30,face="bold"),
          axis.text = element_text(size=28),
          legend.text = element_text(size=28),
          legend.title = element_text(size=30),
          legend.key.height=unit(3,"line"),
          strip.text = element_text(size=28))+
    labs(x="Model",y="Normalised Error" ,colour="Model")+
    scale_y_continuous(limits=c(-1,1))+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  print(P)
  
}
