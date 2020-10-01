rm(list=ls())     # Clean memory
cat("\014")       # Clear Console
graphics.off()    # Close graphs



#install.packages("nleqslv")      # Package to solve a system of non-linear equations.
#install.packages("maxLik")      # Package to optimize the objective function.
library(nleqslv)
library(maxLik)



##########################################################
################ Parameters specification ################
##########################################################

n <- 1000      # sample size.

set.seed(42)      # to make sure that we always generate the same data.

mu <- 0.5      # true value of mu.

sigma <- 2      # standard deviation of the error term

n.sim <- 500      # number of simulations



###########################################
################ Functions ################
###########################################

# function for the maximum likelihood estimation

ML <- function(par, y){      # objective function called ML; "par" an argument of the function and it represents the vector of the unknown parameters that have to be estimated.
  n <- length(y)
  gamma <- par[1]
  theta <- par[2]
  loglike <- rep(NA,n)      # vector to store the observation-specific log-likelihood values.
  for (i in 1:n){
    logl_discr <- 0      # initialize the value for the discrete distribution "part" of the log-likelihood.
    logl_cont <- 0      # initialize the value for the continuous distribution "part" of the log-likelihood.
    if (y[i] == 0) {
      logl_discr <- log(pnorm(-gamma,0,1))       # expression for the discrete distribution "part" of the log-likelihood.
    } else {
      logl_cont <-  -1/2*(log(2*pi) - log(theta**2) +(theta*y[i] - gamma)**2)     # expression for the continuos distribution "part" of the log-likelihood.
    }
    loglike[i] <- logl_discr + logl_cont      # individual (that is, for observation i) log-likelihood.
  }
  return(loglike)      # a vector of observation-specific log-likelihood values must be returned for the optimization function.
}



# function for the method of moments estimation

MM_equations <- function(par, y) {
  moment_eq <- rep(NA,2)      # vector to store the solutions for the moment equations
  mu.hat <- par[1]
  sigma.hat <- par[2]
  # Define the variables from the derivations
  n_1 <- length(y[y>0])
  n <- length(y)
  lambda_h <- dnorm((-mu.hat/sigma.hat))/(1-pnorm(-mu.hat/sigma.hat))
  # insert the two methods of moments equations here:
  moment_eq[1] <- (n_1/n) - pnorm((mu.hat/sigma.hat))
  moment_eq[2] <- 1/n_1 * sum(y) - mu.hat - sigma.hat
  moment_eq
}



#########################################################
################ Monte Carlo simulations ################
#########################################################

# vectors to store the estimated parameters
mu.hat.ML <- rep(NA,n.sim)
sigma.hat.ML <- rep(NA,n.sim)
mu.hat.MM <- rep(NA,n.sim)
sigma.hat.MM <- rep(NA,n.sim)


# Start the simulations

for (s in 1:n.sim){
  
  y_star <- mu + rnorm(n,0,sigma)      # generate the latent variable.
  
  y <- y_star*(y_star > 0)      # observed variable.
  
  init.val <- c(0.0005,0.0002) 
  
  
  # Maximum likelihood estimation:
  
  # Let's optimize using the Newthon-Raphson algorithm:
  optNR <- maxNR(ML, start=init.val, y=y, iterlim=20)      # Maximize the objective function with the Newthon-Raphson method.
  par.hat <- optNR$estimate      # the ML estimates for gamma and theta.
  mu.hat.ML[s] <-   par.hat[1]/par.hat[2]    # the ML estimate for mu.
  sigma.hat.ML[s] <- 1/par.hat[2]       # the ML estimate for sigma.
  
  
  # Method of moments estimation:
  
  result <- nleqslv(init.val, MM_equations, y=y)
  mu.hat.MM[s] <- result$x[1] 
  sigma.hat.MM[s] <- result$x[2] 
  
}


# Compute the mean squared errors for each estimator
MSE.mu.hat.ML= mean((mu.hat.ML-mu)^2)
MSE.sigma.hat.ML= mean((sigma.hat.ML-sigma)^2)
MSE.mu.hat.MM= mean((mu.hat.MM-mu)^2)
MSE.sigma.hat.MM= mean((sigma.hat.MM-sigma)^2)
# outputs for this instance with n=1000:
# > MSE.mu.hat.ML
#[1] 0.004834958
#> MSE.sigma.hat.ML
#[1] 0.003971676
#> MSE.mu.hat.MM
#[1] 0.08212364
#> MSE.sigma.hat.MM
#[1] 1.306659
#Clearly MLE is better, which should make sense, since
# MLE is efficient, so this should give the better estimates (at least quicker)


# print results to screen
MSE.mu.hat.ML
MSE.sigma.hat.ML
MSE.mu.hat.MM
MSE.sigma.hat.MM