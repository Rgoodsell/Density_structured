
# Function to plot observed vs predicted field-level density distributions across models, use with mapply
fDist.check <- function(P,mydata,model){
  
  ## Takes:
  # P: site level transition probabilities to source states [,1:K,] at site [n,,] for all iterations [,,i]
  # mydata: data used in model to generate P
  
  
  # get number of fields & field level ID 
  nField  <- length(unique(mydata$field))
  raneffID <-  as.numeric(mydata$field)
  K <- length(unique(mydata$ytplus))
  
  # set up output objects
  field.diffs <- array(0,dim=c(nField,K)) # df of differences for each field [f,]
  
  
  
  # Loop through fields
  for(f in 1:nField){
    
    obs.dist <- rep(0,5) # observed DS distribution 
    
    
    fsub_preds <-  P[raneffID==f,,] # Select theta estimates by random effect index
    fsub_obs <- mydata[raneffID==f,] # get observed distribution of DS in each field
    
    obsMat <- data.matrix(table(fsub_obs[,2])) # get observed outcome & source dist
    obs <- obsMat/sum(obsMat) # calculate as proportion
    obs.dist[1:length(obs)] <- obs
    
    fsub_mean <-  apply(fsub_preds,2,mean) 
    #fsub_se <- apply(fsub_preds,2,function(x) sqrt(sd(x))/length(x))
    field.diffs[f,] <- fsub_mean -obs.dist  
    
    
  }
  field.diffs  <- as.data.frame.matrix(field.diffs)
  field.diffs$model <- model
  field.diffs <- melt(field.diffs,"model")
  return(field.diffs)
}  


# Plots field-scale distribution error 
plot_field_dist <- function(allDist,levels){
  
  ## Takes:
  # allDist : all errors from all models
  # levels : levels of model variable in alDist, used for reordering
  
  allDist$value <- allDist$value/max(abs(allDist$value))
  distSum <- allDist %>%  group_by(model,variable) %>% summarise(mean = median(value),lower = quantile(value,probs = 0.10),
                                                                 upper = quantile(value,probs=0.90))  
  
  allDist$model <- factor(allDist$model, levels = levels)
  distSum$model <-  factor(distSum$model, levels = levels)
  
  P <- ggplot(distSum,aes(variable,mean))+
    geom_hline(data=distSum,aes(yintercept = 0),size=1,colour="grey",lty=2)+
    geom_jitter(data=allDist,aes(variable,value,colour=model),alpha=0.25,size=4,width=0.1)+
    geom_errorbar(data=distSum,aes(ymax=upper,ymin=lower),size=1,width=0)+
    geom_point(size=5,colour="black") +
    theme_classic()+
    theme(axis.title = element_text(size=30,face="bold"),
          axis.text = element_text(size=28),
          legend.text = element_text(size=28),
          legend.title = element_text(size=30),
          legend.key.height=unit(3,"line"),
          strip.text = element_text(size=28))+
    labs(x="Density State",y="Normalised Error" ,colour="Model")+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    scale_x_discrete(labels= c(1:5))+
    facet_grid(~as.factor(model))
  
  print(P)
  
  
}
