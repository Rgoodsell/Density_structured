##### ---------------- ################################################### ---------------- ######

##### ----------------     Spatial DS modelling   - Covariate functions   ---------------- ######

##### ---------------- ################################################## ---------------- ######



# Gets sum of all surrounding cells from a single x,y coord 
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
  
  
  mx <- mx -1 # rescale back to init ds values
  
  surr.ind <- matrix(surr.ind,ncol=2,nrow=length(surr.ind)/2) # as matrix
  
  
  return(sum(mx[surr.ind])) #  return the sum, -1 to rescale DS states. 
}
### Function to return sums of diagonals and directly adjacent neighbours of a single coord
sum_diags_single <- function(coord,mx){
  
  # coord = x,y coord of cell to be summed
  # mx  = matrix of density states by coordinates in field. Generate using xtabs.
  
  # Get surrounding coords of an individual quadrat
  
  
  # Get diags coords
  diags<- cbind(rep(coord[1] + c(-1,1), times=2), 
                rep(coord[2] + c(-1,1), each=2))
  
  # Get direct coords
  dirs<- cbind(rep(coord[1] + c(-1,0,0,1), times=1), 
               rep(coord[2] + c(0,1,-1,0), each=1))
  
  
  # get an index of out of bounds coords 
  ex<- melt(mx) %>% filter(value==0) %>% select(1:2)  %>% # xtabs creates a rectangular matrix with 0s for empty/unsurveyed cells
    filter(X %in% diags[,1] | X %in% dirs[,1],
           Y %in% diags[,2] | Y %in% dirs[,2])  # Ys
  
  # Get diagonal indeces
  diags.ind  <- diags[!(  diags[,1] < 0| diags[,2] <0| # remove coords if lower than 0 
                            diags[,1] > nrow(mx) | diags[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                            (diags[,1] %in% ex[,1] & diags[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                            (diags[,1]==coord[1] & diags[,2]==coord[2])), ] # remove coords of the cell of interest
  
  
  
  # Get direct indeces
  dir.ind  <- dirs[!(  dirs[,1] < 0| dirs[,2] <0| # remove coords if lower than 0 
                         dirs[,1] > nrow(mx) | dirs[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                         (dirs[,1] %in% ex[,1] & dirs[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                         (dirs[,1]==coord[1] & dirs[,2]==coord[2])), ] # remove coords of the cell of interest
  
  
  mx <- mx -1 # rescale back to init ds values
  
  diag.ind <- matrix(diags.ind,ncol=2,nrow=length(diags.ind)/2) # as matrix
  dir.ind <- matrix(dir.ind,ncol=2,nrow=length(dir.ind)/2) # as matrix
  
  diag.sum <- sum(mx[diag.ind])
  dir.sum <- sum(mx[dir.ind])
  
  out <- data.frame(diag.sum, dir.sum)
  
  return(out) #  return the sum, -1 to rescale DS states. 
}
### Function to return sums of directly adjacent coords, seperating by XY axes, + diagonal coords
sum_adj_single <- function(coord,mx){
  
  # coord = x,y coord of cell to be summed
  # mx  = matrix of density states by coordinates in field. Generate using xtabs.
  
  # Get surrounding coords of an individual quadrat
  
  
  # Get adjacent coords on x axis
  # Get diags corods
  diags<- cbind(rep(coord[1] + c(-1,1), times=2), 
                rep(coord[2] + c(-1,1), each=2))
  
  x_adj<- cbind(coord[1] + c(-1,1), coord[2])
  
  # Get adjacent coords on y axis
  y_adj<- cbind(coord[1], coord[2]+ c(-1,1))
  
  
  # get an index of out of bounds coords 
  ex<- melt(mx) %>% filter(value==0) %>% select(1:2)  %>% # xtabs creates a rectangular matrix with 0s for empty/unsurveyed cells
    filter(X %in% x_adj[,1] | X %in% y_adj[,1] | X %in% diags[,1],
           Y %in% x_adj[,2] | Y %in% y_adj[,2] | Y %in% diags[,2])  # Ys
  
  
  
  # Get x_adjacent indeces
  x.ind  <- x_adj[!(  x_adj[,1] < 0| x_adj[,2] <0| # remove coords if lower than 0 
                        x_adj[,1] > nrow(mx) | x_adj[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                        (x_adj[,1] %in% ex[,1] & x_adj[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                        (x_adj[,1]==coord[1] & x_adj[,2]==coord[2])), ] # remove coords of the cell of interest
  
  
  
  # Get y_adjacent indeces
  y.ind  <- y_adj[!(  y_adj[,1] < 0| y_adj[,2] <0| # remove coords if lower than 0 
                        y_adj[,1] > nrow(mx) | y_adj[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                        (y_adj[,1] %in% ex[,1] & y_adj[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                        (y_adj[,1]==coord[1] & y_adj[,2]==coord[2])), ] # remove coords of the cell of interest
  
  # Get diagonal indeces
  diags.ind  <- diags[!(  diags[,1] < 0| diags[,2] <0| # remove coords if lower than 0 
                            diags[,1] > nrow(mx) | diags[,2] > ncol(mx) |   # remove coords if greater than max X&Y coords
                            (diags[,1] %in% ex[,1] & diags[,2] %in% ex[,2]) | # remove coords of any out of bounds of field
                            (diags[,1]==coord[1] & diags[,2]==coord[2])), ] # remove coords of the cell of interest
  
  
  mx <- mx -1 # rescale back to init ds values
  
  x.ind <- matrix(x.ind,ncol=2,nrow=length(x.ind)/2) # as matrix
  y.ind <- matrix(y.ind,ncol=2,nrow=length(y.ind)/2) # as matrix
  diags.ind <- matrix(diags.ind,ncol=2,nrow=length(diags.ind)/2) # as matrix
  
  
  x.sum <- sum(mx[x.ind])
  y.sum <- sum(mx[y.ind])
  diags.ind <- sum(mx[diags.ind])
  
  out <- data.frame(x.sum, y.sum,diags.ind)
  
  return(out) #  return the sum, -1 to rescale DS states. 
}

# Function to calculate neighbour sums for all coords #
calc_covars <- function(mydata,covar_fun){
  # mydata = DS data 
  # function = numeric index indicating which set of spatial covaraites to be calculated"
  # 1 = isometric neighbour sums
  # 2 = non-isometric - diagonal & direct neighbour sums
  
  ### anon fucntions:
  mx_fun <- function(x){ xtabs(Year1+1~X+Y,x)} # returns matrix of DS states, +1 to fill unsurveyed cells with 0s
  coords_fun <- function(mx){ melt(mx) %>% filter(value>0) %>% as.matrix(.)} # returns matrix of XY coords
  ###
  
  # Split to list by field.year for lapply
  field_list <- split(mydata, mydata$field.year )
  
  mx <- lapply(field_list,mx_fun) # Get matrix
  fx <- lapply(mx,coords_fun) # Get XY coords
  
  # get neighbour sums by field
  
  
  covar_fun <- switch(covar_fun,
                      "1" = sum_neighbour_single,
                      "2" = sum_diags_single,
                      "3" = sum_adj_single)
  
  
  
  ### Need to find a way to retain covar IDs in apply function
  # Or rename after apply (minimizing faff when different covars are being used)
  
  for(i in 1:length(fx)){
    
    xy_split <- split(fx[[i]],seq(nrow(fx[[i]])))
    
    sums <- lapply(xy_split,covar_fun,mx[[i]]) %>% do.call(rbind,.) # apply over xy coords
    
    fx[[i]] <- data.frame(X = fx[[i]][,1], # Xs
                          Y = fx[[i]][,2], # Ys
                          field.year = names(fx)[i]) # field ID
    
    fx[[i]] <- cbind(fx[[i]],sums)
  }
  
  
  fx <- do.call(rbind,fx) # combine
  sum_cov <- left_join(mydata,fx) # Ljoin with original data
  
  return(sum_cov)
  
}



##### ---------------- ############################## ---------------- ######

  ##### ----------------     Data checks    ---------------- ######

##### ---------------- ##############################---------------- ######


### Set of functions to check XY coords rescaled to 1,1 ### 

## returns cases where real xys mismatch desired integer sequence  ## ##
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

### Function that returns fields that change shape between years ## 
check_year_coords <- function(mydata){
  
  # mydata = unscaled data
  
  ## anon functions:
  check_fun <- function(x){
    
    xdiffs <- x %>%  droplevels() %>% pull(X) %>%  split(.,x$Survey)  # get xcoords
    ydiffs <-  x %>%  droplevels() %>% pull(Y) %>% split(.,x$Survey) # ycoords
    
    xs <-  sapply(xdiffs, identical,xdiffs[[1]]) # do they match
    ys <- sapply(ydiffs, identical,ydiffs[[1]]) #
    
    
    return(c(any(xs==F),any(ys==F)))  # return T if no match
    
  }
  ############
  
  # Transform density into a single variable, add in variable for survey year #
  Year1 <- mydata %>% select(-Year2,Crop2) %>% 
    mutate(Density = Year1,Crop=Crop1) %>%
    select(X,Y,Density,Crop,Survey,field)
  
  Year2 <- mydata %>% select(-Year1, -Crop1) %>% 
    mutate(Survey = Survey+1, Density=Year2,Crop=Crop2) %>%
    select(X,Y,Density,Crop,Survey,field)
  
  # combine
  mapping_data <- rbind(Year1,Year2) %>% distinct()
  
  
  # split for lapply
  split <- mapping_data %>% droplevels() %>% split(.,mapping_data$field) 
  
  
# retrieve field names that change shape 
  all <- lapply(split,check_fun) %>% do.call(rbind,.) %>%
    as.data.frame() %>%
    add_rownames() %>% 
    filter(V1 == T, V2 ==T) %>% 
    droplevels() %>% pull(rowname)
  
  return(all)
  
}


##### ---------------- ############################## ---------------- ######

##### ----------------     Plotting Functions   ---------------- ######

##### ---------------- ##############################---------------- ######



## Function to plot a subset of fields over each survey year -  maps changes in bg density ##

map_year_density <- function(rawData,subset){
  
  # rawData = data with rescaled XYs 
  # subset = character vector with names of fields to be mapped

# Transform density into a single variable, add in variable for survey year #
Year1 <- rawData %>% select(-Year2,Crop2) %>% 
  mutate(Density = Year1,Crop=Crop1) %>%
  select(X,Y,Density,Crop,Survey,field)

Year2 <- rawData %>% select(-Year1, -Crop1) %>% 
  mutate(Survey = Survey+1, Density=Year2,Crop=Crop2) %>%
  select(X,Y,Density,Crop,Survey,field)

# combine
mapping_data <- rbind(Year1,Year2) %>% distinct() # remove duplicates - created by splitting data set then +1 to survey (duplicates middle year)


test_dat <- mapping_data %>% filter(field %in% subset)

# Set colours
cols <- c("0" = "grey", "1"="yellow", "2"="orange","3"="red", "4"="dark red")

#plot
map <- ggplot(test_dat,aes(factor(X),factor(Y)))+
  geom_tile(colour="white",aes(fill=factor(Density)))+
  scale_fill_manual(values = cols)+
  facet_grid(Survey ~ field)+
  coord_equal()

print(map)

}