##### ---------------- ################################################### ---------------- ######

##### ----------------     Spatial DS modelling   - Covariate functions   ---------------- ######

##### ---------------- ################################################## ---------------- ######



# Gets sum of surrounding cells from a single x,y coord 
sum_neighbour_single <- function(coord,mx){
  
  # coord = x,y coord of cell to be summed
  # mx  = matrix of density states by coordinates in field. Generate using xtabs.
  
  # Get surrounding coords of an individual quadrat
  surr <- cbind(rep(coord[1] + -1:1, times=3), 
                rep(coord[2] + -1:1, each=3))
  
  
  # get an index of out of bounds coords 
  ex<- melt(mx) %>% filter(value==0) %>% select(1:2)  %>% # xtabs creates a rectangular matrix with 0s for empty/unsurveyed cells
    filter(X %in% surr[,1], # Remove unsurveyed Xs
           Y %in% surr[,2]) # Ys
  
  # if ex is NULL surr.ind breaks 
  surr.ind  <- surr[!(  surr[,1] < 0| surr[,2] <0| # remove coords if lower than 0 
                          surr[,1] > nrow(mx) | surr[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                          (surr[,1] %in% ex[,1] & surr[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                          (surr[,1]==coord[1] & surr[,2]==coord[2])), ] # remove coords of the cell of interest
  
  
  return(sum(mx[surr.ind]-1)) #  return the sum, -1 to rescale DS states. 
}

# Function to calculate neighbour sums for all coords #
calc_neighbour_sums <- function(mydata){
  # mydata = DS date data 
  
  ### anon fucntions:
  mx_fun <- function(x){ xtabs(Year1+1~X+Y,x)} # returns matrix of DS states, +1 to fill unsurveyed cells with 0s
  coords_fun <- function(mx){ melt(mx) %>% filter(value>0) %>% as.matrix(.)} # returns matrix of XY coords
  ###
  
  # Split to list by field.year for lapply
  field_list <- split(myData, myData$field.year )
  
  mx <- lapply(field_list,mx_fun) # Get matrix
  fx <- lapply(mx,coords_fun)     # Get XY coords
  
  # get neighbour sums by field
  for(i in 1:length(fx)){
    
    sums <- apply(fx[[i]],1,sum_neighbour_single,mx[[i]]) # apply over xy coords
    
    fx[[i]] <- data.frame(X = fx[[i]][,1], # Xs
                          Y = fx[[i]][,2], # Ys
                          Sum = sums, # Sums
                          field.year = names(fx)[i]) # field ID
  }
  
  
  fx <- do.call(rbind,fx) # combine
  sum_cov <- left_join(myData,fx) # Ljoin with original data
  
  return(sum_cov)
  
}


##### ---------------- ############################## ---------------- ######

##### ----------------     Rescaled  data checks    ---------------- ######

##### ---------------- ##############################---------------- ######


### Set of functions to check rescaled data ### 

## returns cases where real xys mismatch desired sequence  ## ##
check_diff <- function(mydata){
  
  
  # Mydata = data with XY coords rescaled to start at 1,1
  mydata <- mydata %>% droplevels() 
  
  ## Anon functions ##
  check_fun <- function(x){
    
    # check if difference 
    a <- unique(x$X) # real X
    b <- seq(min(x$X), max(x$X,by=1)) # what X should be
    c <- unique(x$Y) # real Y
    d <- seq(min(x$Y), max(x$Y,by=1)) # what Y should be
    xs<- setdiff(a,b) # X diffs
    ys <- setdiff(c,d) # Y diffs
    
    return(cbind(xs,ys))
    
    
  } # gets difference between real coords & what we want (i.e. integer sequence, without any gaps)
  ######################
  
  # Split into list for lapply
  split <- split(mydata, mydata$field.year)   
  
  

  c<- lapply(split, check_fun) %>%  do.call(rbind,.) # get diffs
  
  if(any(c>0)){ print(c[which(c>0),])} else{print(noquote("All coords start at 1,1"))} # any diffs >1, return rows
  
}

## Function to check if all start at x1,y1  - returns cases where xy != 1,1 ##
check_start <- function(mydata){
  
  # Mydata = data with XY coords rescaled to start at 1,1
  mydata <- mydata %>% droplevels()
  
  #split for lapply
  split <- split(mydata, mydata$field.year)  
  
  
  a<- lapply(split, function(x) min(x$X)) %>%  do.call(rbind,.)# min xs
  b<- lapply(split, function(x) min(x$Y)) %>%  do.call(rbind,.)# min ys
  c<- cbind(a,b)# combine
  
  
  if(any(c!=1)){ print(c[which(c!=1),])} else{print(noquote("All start a X1"))} # if any coords != 1,1 return rows
  
  
}

