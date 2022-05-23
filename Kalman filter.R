
##### Section 3.1 Kalman Filter


### Parameters

n <- 100 # sample size
H <- 1   # Data covariate

true.R <- 0.1 # Data Variance / Measurement error
true.Q <- 0.5 # Parameter Autoregressive Variance
true.M <- 0.7 # AR-1 Term

R <- 0.1
Q <- 0.5
M <- 0.7



### Generating Simulated Data

set.seed(444)
nu <- rnorm(n, 0, true.Q)
x <- rep(NA, n+1)
x[1] <- rnorm(1)

for (t in 2:(n+1)) {
  x[t] <- (true.M)*x[t-1] + nu[t-1]
}

x <- x[-1]

epsilon <- rnorm(n, 0, true.R)
y <- H * x + epsilon

# remove <- c(40:43,80:83)
# remove <- 30:43

# y[remove] <- NA






### Plotting

# plot(x, col = 'orange', type = 'l', ylim = c(min(x)-0.1, max(x)+0.1), xlab = 't')
# points(y, col = 'red', cex = 0.8, pch = 16)
# legend("topright", legend = c("Truth", "Data"),
#        lty = c(1, NA),
#        pch = c(NA, 16),
#        col = c("orange", "red"))







### Filtering

x0 <- 0
P0 <- 10

len.dat <- length(y)

filter.x <- rep(NA, len.dat+1)
filter.Pt <- rep(NA, len.dat+1)

filter.x[1] <- x0
filter.Pt[1] <- P0

forcast.x <- rep(NA, len.dat+1)
forcast.Pt <- rep(NA, len.dat+1)


update.forcast <- function(prev.filter.x, prev.filter.P) {
  forcast.x <- M * prev.filter.x
  forcast.P <- Q + M^2 * prev.filter.P
  return(list(x = forcast.x, P = forcast.P))
}

update.filter <- function(forcast.x, forcast.P, y) {
  K <- forcast.P * H * (H*forcast.P*H + R)^(-1)
  filter.x <- forcast.x + K * (y - H * forcast.x)
  filter.P <- (1-K*H)*forcast.P
  return(list(x = filter.x, P = filter.P))
}

for (t in 2:(len.dat+1)) {
  forcast <- update.forcast(filter.x[t-1], filter.Pt[t-1])
  forcast.x[t] <- forcast$x
  forcast.Pt[t] <- forcast$P
  
  if (is.na(y[t-1])) {
    filter.x[t] <- forcast.x[t]
    filter.Pt[t] <- forcast.Pt[t]
  } else {
    filter <- update.filter(forcast.x[t], forcast.Pt[t], y[t-1])
    filter.x[t] <- filter$x
    filter.Pt[t] <- filter$P
  }
}

# Removing the extra element in front.
forcast.x <- forcast.x[-1]
forcast.Pt <- forcast.Pt[-1]
filter.x <- filter.x[-1]
filter.Pt <- filter.Pt[-1]





### Plotting the filter

plot(x, col = 'orange', type = 'l', ylim = c(min(x)-0.1, max(x)+0.1), xlab = 't')
points(y, col = 'red', cex = 0.8, pch = 16)
lines(filter.x, col = "blue2", lty = 2)
legend("topright", legend = c("Truth", "Data", "Filter"),
       lty = c(1, NA, 2),
       pch = c(NA, 16, NA),
       col = c("orange", "red", "blue2"))



#plot(filter.Pt[2:len.dat], ylim = c(0, max(filter.Pt[2:len.dat])+0.1))










