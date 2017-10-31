##### Ordered category Logit in Stan #####

rm(list=ls())

setwd("/Users/robgoodsell/Google Drive/PhD/R files/DS_model_weeds/Code/Models/STAN/")


library(rstan)
library(dplyr)
library(readr)


###################

##  Data Setup  ##

###################

# Read in data, dop nas & unused levels
rawData <- read.csv("rotation_2007_2008_2009_2010.csv") %>% na.omit(droplevels(.))

str(rawData)

mydata <- rawData %>% 
  filter(Crop2=='wheat', Crop1 =='wheat') %>% 
  mutate(yt = as.factor(Year1+1), ytplus = Year2+1) %>% 
  select(yt, ytplus, field) %>% 
  droplevels(.)



y = mydata$ytplus # Source state (all our data = 1 in this case)
x <- model.matrix(~yt -1, data = mydata)
N = length(mydata$ytplus) # Number of observations
D = dim(x)[2] # Number of explanatory variables 
Fi = length(unique(mydata$field))
field = as.numeric(droplevels(mydata$field))



fixed_code <- read_file("ordered_fixed.stan")




#### Run fixed effects model & Save ####

ordered.logitMod <- stan(model_code =  fixed_code, 
                         data=list(N=N,K=5,D=D,x=x, y = y), warmup = 1000, iter=2000, 
                         chains = 2, cores=2, control = list(adapt_delta = 0.99))



saveRDS(ordered.logitMod,file="fixed_wheatTest")






