
##### Section 5.3.1 Simple Gibbs Sampling Example
### Except I am doing Metropolis-Hastings instead of Gibbs.

library(invgamma)
library(truncnorm)

### Parameters

n <- 100 # sample size
H <- 1   # Data covariate

true.R <- 0.1 # Data Variance / Measurement Error
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

#remove <- c(40:43,80:83)
# remove <- 30:43

#y[remove] <- NA




### Metropolis Hastings


MH <- function(N.mc) {
  # target.d() is target density (up to normalizing constant) function
  # prop.d() is proposal density function
  # ... are arguments to pass to prop.d() and sample.prop()
  
  TF <- c()
  
  
  x <- matrix(NA, ncol = N.mc, nrow = n+1)
  M <- rep(NA, N.mc)
  sigma2eta <- rep(NA, N.mc)
  
  set.seed(444)
  x[,1] <- rnorm(n+1) # Initialization
  M[1] <- 1
  sigma2eta[1] <- 1
  
  for (i in 2:N.mc) {
    #y <- sample.prop(x[i-1], ...)
    #rho <- target.d(y) * prop.d(x[i-1], y, ...) / (target.d(x[i-1]) * prop.d(y, x[i-1], ...))
    #AR <- min(rho,1)
    #TF[i-1] <- sample(c(T, F), 1, prob = c(AR, 1-AR))
    #x[i] <- (TF[i-1] * y) + ((!TF[i-1]) * x[i-1])
    
    prop <- rand.walk.prop(x[,i-1], M[i-1], sigma2eta[i-1])
    #rho <- eval.density(y, prop$x, prop$M, prop$sigma2eta) * rand.walk.dens(x[,i-1], M[i-1], sigma2eta[i-1], prop$x, prop$M, prop$sigma2eta) /
    #  (eval.density(y, x[,i-1], M[i-1], sigma2eta[i-1]) * rand.walk.dens(prop$x, prop$M, prop$sigma2eta, x[,i-1], M[i-1], sigma2eta[i-1]))
    
    numerator <- eval.density(y, prop$x, prop$M, prop$sigma2eta)
    denominator <- eval.density(y, x[,i-1], M[i-1], sigma2eta[i-1])
    
    if (numerator < 1e-10) {
      rho <- 0
    } else {
      rho <- numertor / denominator
    }
    
    #AR <- min(rho,1)
    TF[i-1] <- runif(1) < rho #sample(c(T, F), 1, prob = c(AR, 1-AR))
    
    if (TF[i-1]) {
      x[,i] <- prop$x
      M[i] <- prop$M
      sigma2eta[i] <- prop$sigma2eta
    } else {
      x[,i] <- x[,i-1]
      M[i] <- M[i-1]
      sigma2eta[i] <- sigma2eta[i-1]
    }
  }
  return(list(x.path = x, M.path = M, sigma2eta.path = sigma2eta, accept.rate = mean(TF)))
}



eval.density <- function(y, x, M, sigma2eta) { # Target density
  # y is vector of data
  # x is vector
  # Else are numbers
  
  dens <- prod(dnorm(y[1:n], x[2:(n+1)], R) * dnorm(x = x[2:(n+1)], mean =  M*x[1:n], sd = sqrt(sigma2eta)) * dnorm(x[1]) * dunif(M, -1, 1) * dinvgamma(sigma2eta, 2, scale = 1))
  return(dens)
}


v <- c(0.5, 0.1, 0.1) # step sizes for random walk

rand.walk.prop <- function(x, M, sigma2eta) { # Proposal
  n <- length(x)
  x.new <- rnorm(n+1, mean = x, sd = v[1])
  M.new <- rnorm(1, mean = M, sd = v[2])
  sigma2eta.new <- rtruncnorm(1, a = 0, mean = sigma2eta, sd = v[3])
  return(list(x = x.new, M = M.new, sigma2eta = sigma2eta.new))
}


rand.walk.dens <- function(x.new, M.new, sigma2eta.new, x, M, sigma2eta) { # Proposal density
  dens <- prod(dnorm(x.new, x, v[1])) * dnorm(M.new, M, v[2]) * dtruncnorm(sigma2eta.new, a = 0, mean = sigma2eta, sd = v[3])
  return(dens)
}


mc <- MH(1000)


