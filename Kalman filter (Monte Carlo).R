
### Parameters

n <- 100 # sample size
R <- 0.1 # Data Variance
Q <- 0.5 # Parameter Autoregressive Variance
M <- 0.7 # AR-1 Term
H <- 1   # Data covariate

true.R <- R #0.1
true.Q <- Q #0.5
true.M <- M #0.9



### Generating Simulated Data

set.seed(123)
nu <- rnorm(n, 0, true.Q)
x <- rep(NA, n)
x[1] <- rnorm(1, 0 , 1)

for (t in 2:n) {
  x[t] <- (true.M)*x[t-1] + nu[t]
}

epsilon <- rnorm(n, 0, true.R)
y <- H * x + epsilon

#remove <- c(40:43,80:83)
remove <- 30:43

#y[remove] <- NA


