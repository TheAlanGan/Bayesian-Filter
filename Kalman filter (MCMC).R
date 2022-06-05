
##### Section 5.3.1 Simple Gibbs Sampling Example
### Except I am doing Metropolis-Hastings instead of Gibbs.

library(invgamma)
library(truncnorm)
library(MASS)

### Parameters

n <- 20 # sample size
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
x <- rep(NA, n+1)
x[1] <- rnorm(1)
x0.true <- x[1]

for (t in 2:(n+1)) {
  x[t] <- (true.M)*x[t-1] + nu[t-1]
}

x <- x[-1]
x.true <- x

epsilon <- rnorm(n, 0, true.R)
y <- H * x + epsilon

#remove <- c(40:43,80:83)
# remove <- 30:43

remove <- 10:15
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
  omega <- rep(NA, N.mc)
  
  #set.seed(444)
  x[1,1] <- x0.true
  x[2:(n+1),1] <- x.true #rnorm(n+1, 0, 0.1) # Initialization
  M[1] <- 0.7
  sigma2eta[1] <- 0.5
  
  for (i in 2:N.mc) {

    # Proposal step
    prop <- rand.walk.prop(x[,i-1], M[i-1], sigma2eta[i-1])  # Normal random walk
    # prop <- rand.walk.unif.prop(x[,i-1], M[i-1], sigma2eta[i-1])  # Uniform random walk
    
    # Acceptance probability step
    
    #rho <- eval.density(y, prop$x, prop$M, prop$sigma2eta) * rand.walk.dens(x[,i-1], M[i-1], sigma2eta[i-1], prop$x, prop$M, prop$sigma2eta) /
    #  (eval.density(y, x[,i-1], M[i-1], sigma2eta[i-1]) * rand.walk.dens(prop$x, prop$M, prop$sigma2eta, x[,i-1], M[i-1], sigma2eta[i-1]))
    
    
    ### Note: 'y' here is the observed data.
    # numerator <- eval.density(y, prop$x, prop$M, prop$sigma2eta)
    # denominator <- eval.density(y, x[,i-1], M[i-1], sigma2eta[i-1])
    # 
    # if (numerator < 1e-20) {
    #   rho <- 0
    # } else {
    #   rho <- numerator / denominator
    # }
    # 
    # TF[i-1] <- runif(1) < rho
    #### This is NOT working!!!
    
    
    
    omega[i] <- ratio.of.target.dens(y, prop$x, prop$M, prop$sigma2eta, x[,i-1], M[i-1], sigma2eta[i-1]) *
      rand.walk.dens(x[,i-1], M[i-1], sigma2eta[i-1], prop$x, prop$M, prop$sigma2eta) /
      rand.walk.dens(prop$x, prop$M, prop$sigma2eta, x[,i-1], M[i-1], sigma2eta[i-1])
    TF[i-1] <- runif(1) < omega[i]
    #### This IS working!!!
    
    
    
    # omega[i] <- log.ratio.of.target.dens(y, prop$x, prop$M, prop$sigma2eta, x[,i-1], M[i-1], sigma2eta[i-1])
    # TF[i-1] <- log(runif(1)) < omega[i] # for log.ratio...
    
    
    
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
  return(list(x.path = x, M.path = M, sigma2eta.path = sigma2eta, accept.rate = mean(TF), omega = omega))
}


# Target density
eval.density <- function(y, x, M, sigma2eta) {
  # y is vector of data
  # x is vector
  # Else are numbers
  n <- length(x) - 1
  
  dens <- prod(dnorm(y[1:n], x[2:(n+1)], R) * dnorm(x = x[2:(n+1)], mean =  M*x[1:n], sd = sqrt(sigma2eta)) * dnorm(x[1]) * dunif(M, -1, 1) * dinvgamma(sigma2eta, 2, scale = 1))
  return(dens)
}

# Ratio of target densities
ratio.of.target.dens <- function(y, prop.x, prop.M, prop.sigma2eta, x, M, sigma2eta) {
  # x, M, and sigma2eta are the CURRENT values in the MCMC chain
  # prop.x, prop.M, prop.sigma2eta are the PROPOSED values
  # y is the observed data
  n <- length(x)-1

  omega <- (prop.sigma2eta / sigma2eta)^(-n/2 - 3) *
    exp(1/sigma2eta - 1/prop.sigma2eta) *
    exp(1/2 * (x[1]^2 - prop.x[1]^2)) * 
    exp(1/(2*R) * ( sum((y-x[-1])^2, na.rm = T) - sum((y-prop.x[-1])^2, na.rm = T) )) *
    exp(1/2 * (1/sigma2eta * sum((x[-1]-M*x[-(n+1)])^2) - 1/prop.sigma2eta * sum((prop.x[-1]-prop.M*prop.x[-(n+1)])^2)) )
  
  return(omega)
}

log.ratio.of.target.dens <- function(y, prop.x, prop.M, prop.sigma2eta, x, M, sigma2eta) {
  # x, M, and sigma2eta are the CURRENT values in the MCMC chain
  # prop.x, prop.M, prop.sigma2eta are the PROPOSED values
  # y is the observed data
  n <- length(x)-1
  
  log.omega <- (-n/2 - 3)*log(prop.sigma2eta / sigma2eta) +
    (1/sigma2eta - 1/prop.sigma2eta) +
    (1/2 * (x[1]^2 - prop.x[1]^2)) + 
    (1/(2*R) * ( sum((y-x[-1])^2, na.rm = T) - sum((y-prop.x[-1])^2, na.rm = T) )) +
    (1/2 * (1/sigma2eta * sum((x[-1]-M*x[-(n+1)])^2) - 1/prop.sigma2eta * sum((prop.x[-1]-prop.M*prop.x[-(n+1)])^2)) )
  
  return(log.omega)
}




## Normal random walk
v <- c(1/20, 0.2, 0.2) # step sizes for random walk

rand.walk.prop <- function(x, M, sigma2eta) { # Proposal
  n <- length(x)-1
  
  ### old
  # x.new <- rnorm(n+1, mean = x, sd = v[1])
  # M.new <- rnorm(1, mean = M, sd = v[2])
  # sigma2eta.new <- rtruncnorm(1, a = 0, mean = sigma2eta, sd = v[3])
  
  
  ### new
  cov.mat <- proposal.cov.structure(M, sigma2eta, n) * v[1]
  x.new <- mvrnorm(n=1, mu = x, Sigma = cov.mat)
  
  #M.new <- rnorm(1, mean = M, sd = v[1])
  M.new <- rtruncnorm(1, a = -1, b = 1, mean = M, sd = v[2])
  sigma2eta.new <- rtruncnorm(1, a = 0, mean = sigma2eta, sd = v[3])
  
  return(list(x = x.new, M = M.new, sigma2eta = sigma2eta.new))
}


rand.walk.dens <- function(x.new, M.new, sigma2eta.new, x, M, sigma2eta) { # Proposal density
  #dens <- prod(dnorm(x.new, x, v[1])) * dnorm(M.new, M, v[2]) * dtruncnorm(sigma2eta.new, a = 0, mean = sigma2eta, sd = v[3])
  dens <- dtruncnorm(M.new, a = -1, b = 1, mean = M, sd = v[2]) * dtruncnorm(sigma2eta.new, a = 0, mean = sigma2eta, sd = v[3])
  return(dens)
}





proposal.cov.structure <- function(M, Q, n) {
  Sigma <- matrix(0, nrow = n+1, ncol = n+1)
  
  M.powers <- M^seq(0, 2*n, by = 2)
  vars <- c(1)
  for (i in 1:n) {
    vars[i+1] <- Q * sum(M.powers[0:(i-1)]) + M.powers[i]
  }
  
  diag(Sigma) <- vars
  
  for (j in 1:n) {
    for (i in (j+1):(n+1)) {
      Sigma[i,j] <- M^(i-j) * Sigma[j,j]
      Sigma[j,i] <- Sigma[i,j]
    }
  }
  
  if (!is.positive.definite(Sigma)) {
    stop("Error!! Sigma not pos def!! wtf")
  }
  
  return(Sigma)
}










### Running the MCMC


N.mc <- 200000
burn <- 50000

mc <- MH(N.mc)
mc$accept.rate

mcmc.x <- mc$x.path[,burn:N.mc]
mcmc.M <- mc$M.path[burn:N.mc]
mcmc.sigma2eta <- mc$sigma2eta.path[burn:N.mc]


# plot(mcmc.x[5,], type = 'l')
# hist(mcmc.x[5,], breaks = 50); abline(v=x[5], col = 'red')
# acf(mcmc.x[5,], lag.max = 100)

# plot(mcmc.M, type = 'l')
# hist(mcmc.M, breaks = 30); abline(v = true.M, col = 'red')
# acf(mcmc.M, lag.max = 100)

# plot(mcmc.sigma2eta, type = 'l')
# hist(mcmc.sigma2eta, breaks = 50); abline(v = true.Q, col = 'red')
# acf(mcmc.sigma2eta, lag.max = 100)

filter.x.mcmc <- apply(mcmc.x, 1, mean)
filter.x.2p5 <- apply(mcmc.x, 1, quantile, probs = c(0.025))
filter.x.97p5 <- apply(mcmc.x, 1, quantile, probs = c(0.975))

plot(x.true, col = 'orange', type = 'l', ylim = c(min(x)-0.3, max(x)+0.3), xlab = 't')
points(y, col = 'red', cex = 0.8, pch = 16)
lines(filter.x.mcmc[-1], col = "blue2", lty = 2)
lines(filter.x.2p5[-1], col = "lightblue", lty = 2)
lines(filter.x.97p5[-1], col = "lightblue", lty = 2)
legend("topright", legend = c("Truth", "Data", "Filter"),
       lty = c(1, NA, 2),
       pch = c(NA, 16, NA),
       col = c("orange", "red", "blue2"))

mc$accept.rate

