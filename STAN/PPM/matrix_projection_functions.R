
#-------------------------------------------- ########## Matrix Projection Functions ############# ----------------------------------------------#


#-------------------------------------------- #################################################### ----------------------------------------------#

#### Function to project a single time step with a single matrix ####
project_matrix <- function(matrix,Nt){
  
  # Matrix = matrix to be projected 
  # Nt = matrix of initial ds distributions by row
  
  Nt1 <- apply(Nt,1, function(x) matrix %*% x  )
  
  
  
  return(t(Nt1))
  
  
}


#### Function to project a pair of matrices for a series of time steps ####
project_ts <- function(matrix_pair, steps, Nt){
  
  # Matrix_pair = list of two matrices to be rotated
  # Steps = number of times for rotation to happen
  # Nt = initial DS distribution
  
  names(matrix_pair) <- NULL
  
  # Markov model
  
  ts <- Nt
  
  for(i in 1:steps){
    
    Nt1 <-  matrix_pair[[1]] %*% matrix(Nt) ; ts <- cbind(ts,Nt1)
    Nt <-   matrix_pair[[2]] %*% matrix(Nt1) ; ts <- cbind(ts,Nt)
    
    
    
  }
  return(ts)
  
}


##### Function to project all permutations of fields-level matrices in a single rotation for one full rotation #########
project_permutations_single <- function(rotation_list,Nt){
  
  # Rotation_list = nested list of pairs of rotations
  # Nt = matrix of initial ds distributions by row
  
  library(plyr)
  
  # split arrays of matrices into lists
  step1 <- alply(rotation_list[[1]][,,,1],.margins = 3) # for first half
  step2 <- alply(rotation_list[[2]][,,,1],.margins=3) # second half of rotatioon
  
  # Markov chain model
  iter1 <- lapply(step1,project_matrix,Nt) %>% do.call(rbind,.)  # Project the model for the first step
  iter2 <- lapply(step2,project_matrix,Nt=iter1) %>% do.call(rbind,.) %>% as.data.frame(.)   # & Second step
  
  
  iter2$rotation <- factor(paste(names(rotation_list),collapse=""))  # paste in rotation names
  
  return(iter2)
  
  detach(package:plyr)
  
}


######## Funcion to project all permutations within all rotations + get summary statistics #######
permute_all_single <- function(all_rotations,Nt,states){
  
  #all_rotations  = nested lists of all rotations to be used
  # Nt = matrix of initial DS distributions by row
  
  
  # Project all matrices
  all_permutations   <- lapply(all_rotations,project_permutations,Nt) %>% do.call(rbind,.)
  
  
  # Get summary statistics
  mds <- apply(all_permutations[,1:5],1,mean_ds,states) %>% do.call(rbind,.) # Mean density state
  gma <- apply(all_permutations[,1:5],1,opt.gamma.fun, breaks=breaks, init=c(1,1)) %>% do.call(rbind,.) # Integrate a gamma function across distribution
  
  
  pop.dat <- cbind(mds,gma) # Combine sum stats
  pop.dat$rotation <- all_permutations$rotation # add in rotation ID
  
  return(pop.dat)
  
}


##### returns all permutations of net matrices in a rotation ######
permute_net_stable <- function(matrix_pair){
  
  # matrix_pair = nested list of field-level matrices to be combined
  
  # Internal functions #
  calc_net_matrix <-function(step1,step2){
    
    
    
    return(step1 %*% step2)
    
  } # creates a single net matrix
  #######################
  
  
  library(plyr)
  
  # split arrays of matrices into lists
  step1 <- alply(matrix_pair[[1]][,,,1],.margins = 3) # for first half
  step2 <- alply(matrix_pair[[2]][,,,1],.margins=3) # second half of rotatioon
  
  detach(package:plyr)
  
  out_list <- vector("list",length(step1)) # empty list for loop
  
  for(i in 1:length(step2)) {
    
    out_list[[i]]  <- lapply(step1, FUN=calc_net_matrix, step2[[i]]) # apply net_function for each combination
    
  }
  
  # flatten nesting of lists
  L <- do.call(rbind , out_list)
  # coerce to array
  out <- array(unlist(L), dim = c(nrow(L[[1]]), ncol(L[[1]]), length(L)))
  
  return(  out  ) 
  
}




