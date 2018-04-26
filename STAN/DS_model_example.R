#-----------------################################################ ----------------------------------------------################################################ -----------------------------#

#-------------------------------------------------------------------########### Full Example of DS model ########### ---------------------------------------------------------------------------#

#-----------------################################################ ----------------------------------------------################################################ -----------------------------#

# Clear workspace
rm(list=ls()) 


# Setwd
setwd("/Users/robgoodsell/Google Drive/PhD/R_files/DS_model_weeds/Code/Models/STAN/") 


# Load packages 
required_pkgs <- c( "ggplot2", "dplyr" , "rstan" , "readr" ) #
if( !require(pacman) ) install.packages( "pacman" )
pacman::p_load(required_pkgs, character.only = T)



#------------------- ################### -------------------------------------- ################### -------------------------------------- ################### -------------------#

#------------------- ################### -------------------------------------- ##  Data Setup  ##  ------------------------------------- ################### -------------------#

#------------------- ################### ------------------------------------- ################### -------------------------------------- ################### --------------------#


# Read in raw data - only includes fields that have black-grass infestation at some point... 
rawData <- readRDS("rotation_allstates.rds") %>% na.omit(droplevels(.))




### Set up data ####
mydata <- rawData %>%
  mutate(yt =Year1+1, ytplus = (Year2+1)) %>%  # get rid of 0 states
  filter(Crop1 == "wheat", # Select rotation
         Crop2 == "wheat") %>% 
  select(yt, ytplus, field=field.year,Crop1,Crop2) %>% # select density states, crops & 
  droplevels(.) 





## Check data ##
str(mydata)

# yt = density state first year
# ytplus = density state second year
# field = field ID
# Crop 1 & 2 = crops in year 1 & 2





##### Set up Stan Data ########

## Population level ##
y = mydata$ytplus; # Destination state / outcome 
K = 5 # Number of density states
x <- model.matrix(~factor(yt) -1, data = mydata); # Design matrix of source state as primary covariate 
N = length(mydata$ytplus); # Number of observations
D = dim(x)[2]; # Number of explanatory variables 



# Group-level effects for models II - V #
Fi = length(unique(mydata$field)); # Number of fields / group-level effects
field = as.numeric(droplevels(mydata$field)) # indexing variable ; recode to numeric



## Hyper-parameters for covariance matrix - only used for model III ##
# lkj_param = 1; # Set priors for covariance matrix 
# sd_diags =rep(1.25,5)



#------------------- ################### -------------------------------------- ################### -------------------------------------- ################### -------------------#

#------------------- ################### -------------------------------------- ##  Stan model  ##  ------------------------------------ ################### --------------------#

#------------------- ################### -------------------------------------- ################### -------------------------------------- ################### -------------------#

## Select which model to run ## 

# Model files 
list.files("Stan files/",pattern = ".stan$")

# ordered_fixed.stan                = fixed effects model (Model I) 
# ordered_random_intercept.stan     = single group level effect (Model II)
# ordered_multi_cholesky_unif.stan  = group level source state effects (Model III)
# ordered_random_cutpoints.stan     = group level cutpoints (Model IV)
# ordered_random_cutints.stan       = group level cutpoints + single intercept (Model V)



# Decide on model # 
stanCode <- read_file("Stan files/ordered_random_intercept.stan")



# Combine relevant data into a list # 
stanData <-  list(y = y, N = N, K = K,  x = x, D = D, Fi = Fi, field = field)




## Run the model ##
ordered.raneffMod <- stan (  model_code =  stanCode, # Stan code
                                  data  =  stanData, # Stan data
                  
                                warmup  = 1000, # Sampling parameters
                                  iter  = 2000, 
                                chains  = 4,
                                 cores  = 4, # N cores for paralellisation
                          
                               control  = list(adapt_delta=0.99,   # Lowering adapt_delta (similar to step size in MH MCMC) can improve sampling times but reduces efficiency of the algorithim to explore the posterior,
                                                                  # generally needs to be above 0.9 to avoid step-size warnings in these models 
                                              
                        max_treedepth   = 15)) # Tree depth can also be lowered but usually throws errors lower than 15






#------------------- ################### -------------------------------------- ################### -------------------------------------- ################### -------------------#

#------------------- ################### -------------------------------------- ##  Diagnostics  ##  ------------------------------------ ################### --------------------#

#------------------- ################### -------------------------------------- ################### -------------------------------------- ################### -------------------#




# Check model parameters #
ordered.raneffMod@model_pars



# Model summary #
summary(ordered.raneffMod)



# Inspect population level effects
plot(ordered.raneffMod,pars=c("beta","c"))
traceplot(ordered.raneffMod,pars=c("beta","c"))



# Inspect Group level parameters (pars will change depending on model ran)
plot(ordered.raneffMod,pars="gamma")
traceplot(ordered.raneffMod,pars="sigmaGamma") # Check hyperparameters
traceplot(ordered.raneffMod,pars=c("gamma")) # Check subsets if lots of groups


## Check effective sample size and scale reduction factors ##

## Whole model ##
stan_ess(ordered.raneffMod) # ESS
stan_rhat(ordered.raneffMod) #rhat


# Individual parameters #
stan_ess(ordered.raneffMod,pars="beta") # ESS
stan_rhat(ordered.raneffMod,pars="beta") #rhat



## Save model ##
#saveRDS(ordered.raneffMod, "stanMod.rds")


