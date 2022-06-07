
##### Section 5.3.1 Simple Gibbs Sampling Example
### Doing Gibbs using the full conditionals derived in the paper.


library(invgamma) # For inverse gamma
library(truncnorm) # For truncated normal
library(MASS)
library(matrixcalc) # For is.positive.def()

### Parameters

n <- 15 # sample size
H <- 1   # Data covariate

true.R <- 0.1 # Data Variance / Measurement Error
true.Q <- 0.5 # Parameter Autoregressive Variance
true.M <- 0.7 # AR-1 Term

R <- 0.1
Q <- 0.5
M <- 0.7


### Generating Simulated Data

set.seed(123)
nu <- rnorm(n, 0, sd = sqrt(true.Q))
x <- rep(NA, n+1)
x[1] <- rnorm(1)
x0.true <- x[1]

for (t in 2:(n+1)) {
  x[t] <- (true.M)*x[t-1] + nu[t-1]
}

x <- x[-1]
x.true <- x

epsilon <- rnorm(n, 0, sd = sqrt(true.R))
y <- H * x + epsilon

# remove <- c(40:43,80:83)
# remove <- 30:43

# remove <- 10:15
# y[remove] <- NA








Gibbs <- function(N.mc, y) {
  
  n <- length(y)

  x <- matrix(NA, ncol = N.mc, nrow = n+1)
  M <- rep(NA, N.mc)
  sigma2eta <- rep(NA, N.mc)
  omega <- rep(NA, N.mc)
  
  #set.seed(444)
  x[1,1] <- x0.true
  x[2:(n+1),1] <- x.true #rnorm(n+1, 0, 0.1) # Initialization
  M[1] <- 0.7
  sigma2eta[1] <- 0.5
  
  
  for (i in 2:N.mc) {
    x[1, i] <- rnorm(1, mean = (M[i-1]*x[2,i-1]/sigma2eta[i-1])/(M[i-1]^2/sigma2eta[i-1] + 1), sd = sqrt(1/(M[i-1]^2/sigma2eta[i-1] + 1)) )
    
    for (j in 2:n) {
      
      if (!is.na(y[j-1])) {
        x[j, i] <- rnorm(1, mean = (M[i-1]*x[j-1, i]/sigma2eta[i-1] + M[i-1]*x[j+1, i-1]/sigma2eta[i-1] + y[j-1]/R)/(M[i-1]^2/sigma2eta[i-1] + 1/sigma2eta[i-1] + 1/R), sd = sqrt(1/(M[i-1]^2/sigma2eta[i-1] + 1/sigma2eta[i-1] + 1/R)))
      } else {
        x[j, i] <- rnorm(1, mean = (M[i-1]*x[j-1, i]/sigma2eta[i-1] + M[i-1]*x[j+1, i-1]/sigma2eta[i-1])/(M[i-1]^2/sigma2eta[i-1] + 1/sigma2eta[i-1]), sd = sqrt(1/(M[i-1]^2/sigma2eta[i-1] + 1/sigma2eta[i-1])))
      }
    }
    
    if (!is.na(y[n])) {
      x[(n+1), i] <- rnorm(1, mean = (M*x[n, i]/sigma2eta[i-1] + y[n]/R)/(1/sigma2eta[i-1] + 1/R), sd = sqrt(1/(1/sigma2eta[i-1] + 1/R)))
    } else {
      x[(n+1), i] <- rnorm(1, mean = (M*x[n, i]/sigma2eta[i-1])/(1/sigma2eta[i-1]), sd = sqrt(1/(1/sigma2eta[i-1])))
    }
    
    M[i] <- rtruncnorm(a = -1, b = 1, n = 1, mean = sum(x[2:(n+1), i] * x[1:n, i] / sigma2eta[i-1]) / sum(x[1:n, i]^2 / sigma2eta[i-1]), sd = sqrt(1/(sum(x[1:n, i]^2 / sigma2eta[i-1]))))
    
    sigma2eta[i] <- rinvgamma(1, shape = n/2 + 2, scale = 1/(1 + 0.5 * sum((x[2:(n+1), i] - M[i] * x[1:n, i])^2)))
    
    if (i %% round(0.05*N.mc) == 0) { # To monitor progress
      cat(round(i/N.mc, 2)*100, "% done \n", sep = "")
    }
  }

  return(list(x.path = x, M.path = M, sigma2eta.path = sigma2eta))
}





### Running the gibbs


N.mc <- 100000
burn <- 10000

mc.gibbs <- Gibbs(N.mc, y)

gibbs.x <- mc.gibbs$x.path[,burn:N.mc]
gibbs.M <- mc.gibbs$M.path[burn:N.mc]
gibbs.sigma2eta <- mc.gibbs$sigma2eta.path[burn:N.mc]


# plot(gibbs.x[5,], type = 'l')
# hist(gibbs.x[5,], breaks = 50, probability = T); abline(v=x[5], col = 'red'); abline(v = mean(gibbs.x[5,]), col = 'blue', lty = 2); lines(density(gibbs.x[5,], n = 50))
# acf(gibbs.x[5,], lag.max = 100)
# 
# plot(gibbs.M, type = 'l')
# hist(gibbs.M, breaks = 30, probability = T); abline(v = true.M, col = 'red'); abline(v = mean(gibbs.M), col = 'blue', lty = 2); lines(density(gibbs.M, n = 512))
# acf(gibbs.M, lag.max = 100)

# plot(gibbs.sigma2eta, type = 'l')
# hist(gibbs.sigma2eta, breaks = 50, probability = T); abline(v = true.Q, col = 'red'); abline(v = mean(gibbs.sigma2eta), col = 'blue', lty = 2)
# acf(gibbs.sigma2eta, lag.max = 100)

filter.x.gibbs <- apply(gibbs.x, 1, mean)
filter.x.2p5.gibbs <- apply(gibbs.x, 1, quantile, probs = c(0.025))
filter.x.97p5.gibbs <- apply(gibbs.x, 1, quantile, probs = c(0.975))

plot(x.true, col = 'orange', type = 'l', ylim = c(min(x,y,na.rm=T)-0.3, max(x,y,na.rm=T)+0.3), xlab = 't', main = "Gibbs")
points(y, col = 'red', cex = 0.8, pch = 16)
lines(filter.x.gibbs[-1], col = "blue2", lty = 2)
lines(filter.x.2p5.gibbs[-1], col = "lightblue", lty = 2)
lines(filter.x.97p5.gibbs[-1], col = "lightblue", lty = 2)
legend("topright", legend = c("Truth", "Data", "Filter"),
       lty = c(1, NA, 2),
       pch = c(NA, 16, NA),
       col = c("orange", "red", "blue2"))


