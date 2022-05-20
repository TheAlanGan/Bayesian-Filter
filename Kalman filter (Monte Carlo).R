
##### Section 4.3.2 Ensemble Kalman Filter


### Parameters

n <- 100 # sample size
m <- 100 # Number of Monte Carlo samples


H <- 1   # Data covariate
true.R <- 0.1 # Data Variance
true.Q <- 0.5 # Parameter Autoregressive Variance
true.M <- 0.7 # AR-1 Term

R <- 0.1
Q <- 0.5
M <- 0.7


### Generating Simulated Data

set.seed(123)
nu <- rnorm(n, 0, true.Q)
x <- rep(NA, n)
x[1] <- rnorm(1)

for (t in 2:n) {
  x[t] <- (true.M)*x[t-1] + nu[t]
}

epsilon <- rnorm(n, 0, true.R)
y <- H * x + epsilon

remove <- c(40:43,80:83)
# remove <- 30:43

y[remove] <- NA




### Defining functions

forward.evolution <- function(prev.filter.x) {
  # Takes in previous analysis/filter value (from MC) and updates it to return
  # next forecast value.
  
  m <- length(prev.filter.x)
  forcast.x <- M*prev.filter.x + rnorm(m, 0, Q)
  return(forcast.x)
}

get.forcast.Pt <- function(forcast.x) {
  # Gets the sample variance of forcasts.x to get sample P.
  
  forcast.Pt <- var(forcast.x)
  return(forcast.Pt)
}


filter.update <- function(forcast.x, forcast.Pt, y) {
  # Takes in current forcast value and uses Kalman filter update equations to 
  # give sampled filter/analysis value.
  
  m <- length(forcast.x)
  K <- forcast.Pt * H * (H * forcast.Pt * H + R)^(-1)
  
  filter.x <- forcast.x + K * (y + rnorm(m,0,R) - H*forcast.x)
  return(filter.x)
}






### The filtering

filter.x <- matrix(NA, nrow = m, ncol = n+1)
forcast.x <- matrix(NA, nrow = m, ncol = n+1)
P <- rep(NA, n+1)

filter.x[,1] <- rnorm(m, 0, 1)



for (t in 2:(n+1)) {
  forcast.x[,t] <- forward.evolution(filter.x[,(t-1)])
  P[t] <- get.forcast.Pt(forcast.x[,t])
  
  if (is.na(y[t-1])) {
    filter.x[,t] <- forcast.x[,t]
  } else {
    filter.x[,t] <- filter.update(forcast.x[,t], P[t], y[t-1])
  }
  
}


filter.x.mean <- apply(filter.x, 2, mean)
filter.x.mean <- filter.x.mean[-1]
filter.Pt <- apply(filter.x, 2, var)
filter.Pt <- filter.Pt[-1]





### Plotting the filter

plot(x, col = 'orange', type = 'l', ylim = c(min(x)-0.1, max(x)+0.1), xlab = 't')
points(y, col = 'red', cex = 0.8, pch = 16)
lines(filter.x.mean, col = "blue2", lty = 2)
legend("topright", legend = c("Truth", "Data", "Filter"),
       lty = c(1, NA, 2),
       pch = c(NA, 16, NA),
       col = c("orange", "red", "blue2"))



#plot(filter.Pt, type = 'l', ylim = c(0, 1.4))













