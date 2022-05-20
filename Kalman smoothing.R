
source("C:/Users/Alang/Downloads/Summer - 2022 - OC/Kalman filter.R")

### Parameters

# n <- length(y) # sample size
# R <- 0.1 # Data Variance
# Q <- 0.5 # Parameter Autoregressive Variance
# M <- 0.7 # AR-1 Term
# H <- 1 # Data covariate


J <- rep(NA, n)
smooth.x <- rep(NA, n)
smooth.Pt <- rep(NA, n)

smooth.x[n] <- filter.x[n]
smooth.Pt[n] <- filter.Pt[n]


for (t in (n-1):1) {
  J[t] <- filter.Pt[t] * M * forcast.Pt[t+1]^(-1)
  smooth.x[t] <- filter.x[t]  + J[t] * (smooth.x[t+1] - forcast.x[t+1])
  smooth.Pt[t] <- filter.Pt[t] + J[t] * (smooth.Pt[t+1] - forcast.Pt[t+1]) * J[t]
}




### Plotting the smoother

plot(x, col = 'orange', type = 'l', ylim = c(min(x)-0.1, max(x)+0.1), xlab = 't')
points(y, col = 'red', cex = 0.8, pch = 16)
lines(smooth.x, col = "limegreen", lty = 2)
legend("topright", legend = c("Truth", "Data", "Smoother"),
       lty = c(1, NA, 2),
       pch = c(NA, 16, NA),
       col = c("orange", "red", "limegreen"))





### smoother vs. filter

plot(x, col = 'orange', type = 'l', ylim = c(min(x)-0.1, max(x)+0.1), xlab = 't')
points(y, col = 'red', cex = 0.8, pch = 16)
lines(smooth.x, col = "limegreen", lty = 2)
lines(filter.x, col = "blue2", lty = 2)
legend("topright", legend = c("Truth", "Data", "Smoother", "Filter"),
       lty = c(1, NA, 2, 2),
       pch = c(NA, 16, NA, NA),
       col = c("orange", "red", "limegreen", "blue2"))





