rm(list=ls())
setwd("C:/Users/Robert/Documents/PhD/R files/DS_model_weeds/Iceberg_scripts/results")# set path to the distribution folder

library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(scales)
library(gridExtra)
library(magrittr)


####### Transition matrix heatmaps ########



hm<-function(matrix){
mat<-read.table(file=matrix,sep="\t",dec=".")

mat<-round(mat, digits = 4)

mat$States<-c("1","2","3","4","5")

melted<-melt(mat)

ggplot(melted, aes(variable, States)) + geom_tile(aes(fill = value),colour = "white") + 
  
  geom_text(aes(label=value, fontface="bold"), size=6)+
  
  scale_fill_gradient(low = "white",high = "Dark blue") +
  
  scale_x_discrete(name="State in Year 1",
                   breaks=c("V1", "V2", "V3","V4","V5"),
                   labels=c("A", "L", "M","H","VH"))+
  
  scale_y_discrete(name="State in Year 2",
                   breaks=c("1", "2", "3","4","5"),
                   labels=c("A", "L", "M","H","VH"))+

  theme(axis.title.x = element_text(face="bold", colour="#990000", size=20),
        axis.title.y  = element_text(face="bold", colour="#990000", size=20),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))
        

}

naive_hm <- hm("naive_ib.txt")   + ggtitle("Naive") + theme(plot.title=element_text(size=20, face="bold"))
norot_hm <- hm("norotation.txt") + ggtitle("Unrotated") + theme(plot.title=element_text(size=20, face="bold"))
rot_hm   <- hm("rotation.txt")   + ggtitle("Rotated") + theme(plot.title=element_text(size=20, face="bold"))
year1_hm <- hm("year1.txt")      + ggtitle("First Survey Year") + theme(plot.title=element_text(size=20, face="bold"))
year2_hm <- hm("year2.txt")      + ggtitle("Second Survey Year") + theme(plot.title=element_text(size=20, face="bold"))

grid.arrange(naive_hm,norot_hm,rot_hm,year1_hm,year2_hm, nrow=3, ncol=2)

#######  Stable density structure figures #######

#### two equivalent ways ###


### 1 - direct iteration ####

state.distrib<-function(matrix){ 
  
  av.mat<-read.table(file=matrix,sep="\t",dec=".")
  
  av.mat<-round(av.mat, digits = 4)
  
  av.mat<-as.matrix(t(av.mat))
  
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
    if (counter==2000) break
    
  }
  
  # combine vectors into data frame and plot
  f<-as.vector(vec1)
  s<-as.vector(states)
  final<-data.frame(s,f)
  


  ggplot(final,aes(x=factor(s),y=f)) + geom_bar(stat="identity")  + 
    
    scale_y_continuous(limits = c(0, 1), name="Proportion") + 
    
    scale_x_discrete(name="Density State",
                     breaks=c("1", "2", "3","4","5"),
                     labels=c("A", "L", "M","H","VH"))
    
  
  
  
    
  
}


naive_ssd <- state.distrib("naive_ib.txt") 
norot_ssd <- state.distrib("norotation.txt")
rot_ssd   <- state.distrib("rotation.txt")
year1_ssd <- state.distrib("year1.txt")
year2_ssd <- state.distrib("year2.txt")

grid.arrange(naive_ssd,year1_ssd,year2_ssd,norot_ssd,rot_ssd, nrow=3, ncol=2)


##### 2 - solving eigen system ##### 

state.distrib2 <-function(matrix){
  
  
  mat<-read.table(file=matrix,sep="\t",dec=".")
  
  
  es <- eigen( mat )
  
  n2 <- es$vector[,1] / sum(es$vector[,1])
  
  n2 <-as.numeric(n2)
  
  n2
  
  return(n2)
}



ns<-state.distrib2("naive_ib.txt")
nrs<-state.distrib2("norotation.txt")
rs<- state.distrib2("rotation.txt")
y1s<- state.distrib2("year1.txt")
y2s<- state.distrib2("year2.txt")




#### Summing non-zero states for occupancy ###

final<-data.frame(ns,y1s,y2s,nrs,rs)

finalapply <-apply(final[2:5,], MARGIN=2, FUN=sum)

model<- c("1","2","3","4","5")

df<-data.frame(finalapply,model)

ggplot(df, aes(x=model,y=finalapply)) + geom_bar(stat="identity", width=0.7, position=position_dodge(2)) +
  
  scale_x_discrete(name="Model",
                   breaks=c("1", "2", "3","4","5"),
                   labels=c("N", "Y1", "Y2","NR","R"))+
  
  scale_y_continuous(name="Proportion")+
  
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.title.y  = element_text(face="bold",  size=20),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))



 ### damping ratios #### 
#install.packages("popbio")
library(popbio)

  naive<-read.table("naive_ib.txt")
  norot<-read.table("norotation.txt")
  rot<-read.table("rotation.txt")
  year1<-read.table("year1.txt")
  year2<-read.table("year2.txt")
  
  
matrix.list<-list(naive,year1,year2,norot,rot)


damping.ratio<- function(av.matrix){
  eigen.vals<- eigen(av.matrix)
  DR<-Mod(eigen.vals$values)
  DR2 <- DR[2]
  
}

dr<-lapply(matrix.list,t) %>%
lapply(damping.ratio) %>%
unlist()%>%
data.frame()
dr[,2]<-c("1","2","3","4","5")

names(dr)<-c("DR","Model")

ggplot(dr,aes(x=Model,y=DR)) + geom_bar(stat="identity", colour="black", fill="seagreen", size=1.3)+
  

  scale_x_discrete(name="Model",
                   breaks=c("1", "2", "3","4","5"),
                   labels=c("N", "Y1", "Y2","NR","R"))+
  
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.title.y  = element_text(face="bold",  size=20),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))



#### gamma dist mean & variance ####

setwd("C:/Users/Robert/Documents/PhD/R files/DS_model_weeds/Bayesian model")# set path to the distribution folder


source("DS_integration_code.R")


##### get stable state dist and remove absent proportion #####
naive_states<-c( 0.146640624, 0.022180724, 0.013461759, 0.004526887)
no_rotation_states<-c( 0.11939217, 0.02278082, 0.03122278, 0.01298473)
rotation_states<-c( 0.139779896, 0.017298175, 0.007673011, 0.002672391)

# breaks for gamma func ## 

breaks <-log(c(1,160,450,1450,2000))


# run code #
naive<-optGammaFun(naive_states,breaks)
no_rot<-optGammaFun(no_rotation_states,breaks)
rot<-optGammaFun(rotation_states,breaks)

data<-data.frame(c(naive$mean,no_rot$mean,rot$mean),c(naive$var,no_rot$var,rot$var),c("Naive","Unrotated","Rotated"))
names(data)<-c("Mean","Variance","Model")
data[,2]<-sqrt(data$Variance)
data

GI<-ggplot(data,aes(Model,Mean)) + geom_point() + geom_errorbar(aes(x=Model, ymin=Mean-Variance,ymax=Mean+Variance, width=0.5))+
  
  scale_y_continuous(name="Log Average Mean")+
  
  theme(axis.title.x = element_text(face="bold", size=20),
        axis.title.y  = element_text(face="bold",  size=20),
        axis.text.x = element_text(size=17),
        axis.text.y = element_text(size=17))



save.image(file="outData.Rdata")
