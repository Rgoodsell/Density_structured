# ------- Demonstrates numerical integration of DS data ---- 

rm( list = ls() )

# --- Custom function to create binned data ---
robhist <- function( x, breaks) {		
		xmin <- min(x)#computes the minmum value of a vector
		xmax <- max(x)#computes the maximun value of a vector
		step <- (xmax - xmin) / breaks # takes difference and / by the user defined number of breaks
		intervals <- seq(xmin, xmax, by = step) # uses the size of step to make all the points inbetween the min and max
		countfun <- function( i ) length( intersect( which( x > intervals[i]), which( x < intervals[i+1] ) ) ) 
    # intersection gives the numbers that overlap between the data sets
    #gives back a vector that is then used in the sapply as the function in sapply 
		dx <- sapply( c(1: (length( intervals) - 1) ), countfun)
    #dx gives a vector back not a list, using sapply fuinction 
		return( list( counts = dx, breaks = intervals) )	
    #so counts is a vector of the number of weeds in a plot?
    #breaks is the size along the x axies 
	}


# --- Integrates log-normal ----
lnormFun <- function( breaks, lmean, lvar ) {
		
		n <- length(breaks)	
		starts <- breaks[1:(n-1)] 
		ends <-  breaks[ 2:n]  
		probs <- plnorm( ends, lmean, lvar) - plnorm( starts, lmean, lvar )
    #ends= vector of quantiles, lmean=mean log, lvar= sd log 
    #probs = dumlative density function for the difference in the final break and the start break
		return( probs)
	
	}
	
	
# --- Optimizes the lognormal function for data and given breaks --- 	
optlnormFun <- function( ds, breaks, init = c(1, 1) ) {
	
	# --- Function that is optimised, crudely based on least squares --- 
	optFun <- function( pars) {
		 lmean <- pars[1]
		 lvar <-pars[2]
		 probs <- lnormFun( breaks, lmean, lvar )
		 ssq <-  ( ds - probs )^2
		 return( sum( ssq) )
	}
	
	ds <- ds / sum(ds)
	op <- optim( init, optFun)
	n <- length(breaks)
	lmean <- op$par[1]
	lvar <- op$par[2]
	starts <-  breaks[1:(n-1)] 
	ends <-  breaks[ 2:n] 
	
	probs <- plnorm(ends, lmean, lvar) - plnorm(starts, lmean, lvar)
	
	x <- seq(min(breaks), max(breaks), by = (max(breaks) - min(breaks)) / 100)
	fx <-  plnorm(x[2:length(x) ], lmean, lvar) - plnorm(x[1:(length(x)-1)], lmean, lvar) ###have checked and not needed..
	plot( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, probs, type = "l" , log = "x",
		xlab = "State", ylim = c(0,1))
	points( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, ds, col= "red" )
	mean <- exp( lmean + 0.5 * lvar )
	var <- ( exp( lvar ) - 1 ) * exp( 2 * lmean + lvar )
	
	return( list( lmean = lmean, lvar = lvar, mean = mean, var = var, probs = probs, ds = ds) )	
}


# ---- Integrates Gamma Function --- 
GammaFun <- function( breaks, shape, rate ) {
		n <- length(breaks)
		starts <- breaks[1:(n-1)] 
		ends <-  breaks[ 2:n]  
		probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate)
		return( probs )
	}
	

# --- Optimizes the gamma function for data and given breaks --- 	
optGammaFun <- function( ds, breaks, init = c(1, 1) ) {
	
	# --- Function that is optimised, crudely based on least squares --- 
	optFun <- function( pars) {
		 shape <- pars[1]
		 rate <-pars[2]
		 probs <- GammaFun( breaks, shape, rate)
		 ssq <-  ( ds - probs )^2
		 
		 return( sum( ssq) )
	}
	
	ds <- ds / sum(ds)
	op <- optim( init, optFun)
	n <- length(breaks)
	shape <- op$par[1]
	rate <- op$par[2]
	starts <-  breaks[1:(n-1)] 
	ends <-  breaks[ 2:n] 
	
	probs <- pgamma(ends, shape, rate) - pgamma(starts, shape, rate)
	
	x <- seq(min(breaks), max(breaks), by = (max(breaks) - min(breaks)) / 100)
	fx <-  pgamma(x[2:length(x) ], shape, rate) - pgamma(x[1:(length(x)-1)], shape, rate)
	plot( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, probs, type = "l" , log = "x",
		xlab = "State", ylab = "Proportion", ylim = c(0,1))
	points( ( breaks[1:(length(breaks)-1) ]  + breaks[ 2:length(breaks)  ] ) / 2, ds, col= "red" )
	mean <- shape / rate
	var <- shape / rate^2
	
	return( list(shape = shape, rate = rate, mean = mean, var = var, probs = probs, ds = ds) )
	
}





# ------- run some code ------ 


# Generate data with known parameters
set.seed(1)
par(mfrow=c(2,2))

# the SDlog value dictates the shape of the graph 

x <- rlnorm( 10000, 0, sdlog = 0.2) # the mean log of 0 gives a mean of 1 

# Coarsen the data
breaks <- 5
DSdata <- robhist( x, breaks) 
hist(DSdata$counts)

# Fit a model to the coarsened data
optlnormFun( DSdata$counts, DSdata$breaks)
# for the gamma 
optGammaFun(DSdata$counts, DSdata$breaks)

mean_log =optlnormFun( DSdata$counts, DSdata$breaks)$mean
mean_gamma =optGammaFun(DSdata$counts, DSdata$breaks)$mean
diff_mean_log= mean_log-1
diff_mean_gamma= mean_gamma-1

if (diff_mean_log < diff_mean_gamma) { print("log is closer to the mean")}  else { print("gamma is closer to the mean")}

set.seed(1)
x <- rlnorm( 10000, 0, 0.5 )

bks <- 5
DSdata <- robhist( x, bks) 

optGammaFun( DSdata$counts, DSdata$breaks)




set.seed(1)
x <- rlnorm( 10000, 2, 0.2 )

bks <- 5
DSdata <- robhist( x, bks) 
hist(DSdata$counts)

optGammaFun( DSdata$counts, DSdata$breaks)



set.seed(1)
x <- rlnorm( 10000, 2, 0.5 )

bks <- 5
DSdata <- robhist( x, bks) 

optlnormFun( DSdata$counts, DSdata$breaks)



# This one ^^^!^!^!^ is harder to fit to 5 states
set.seed(1)
x <- rlnorm( 10000, 3, 1.5 ) # originally 10000,3,1.5

bks <- 5
DSdata <- robhist( x, bks) 

optlnormFun( DSdata$counts, DSdata$breaks)


# Needs more states...... 


bks <- 100
DSdata <- robhist( x, bks) 
hist(DSdata$counts)

optlnormFun( DSdata$counts, DSdata$breaks)






