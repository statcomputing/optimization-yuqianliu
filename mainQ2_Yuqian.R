x2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)

## lodaing package
library(ggplot2)

## Problem 2.a loglikelihood
Theta_value2 <- seq(-pi, pi, by=0.01)

loglik2 <- function(theta){
  # theta : location
  # x : a vector with samples
  x2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
          2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
  n <- length(x2)
  term <- 0
  for(i in 1:n){
    term <- term + log(1-cos(x2[i]-theta))
  }
  ll <- -2*n*log(pi)+ term
  ll
}

LogLik_plot2 <- sapply(Theta_value2, function(theta) loglik2(theta))
# df <- data.frame(Theta_value, LogLik_plot)
# g <- ggplot(df, aes(Theta_value, LogLik_plot))
# plot log-likelihood
plot(Theta_value2, LogLik_plot2,xlab = "theta",ylab = "log-likelihood",
     main = "Log-likelihood Function")
#---test---#
#likelihood <- function(theta){
#  x2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
#          2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
#  n <- length(x2)
#  L1 <- 1
#  for(i in 1:n){
#    L1 <- L1*(1-cos(x2[i]-theta))/(2*pi)
#  }
#  L1
#}
#Lik_plot2 <- sapply(Theta_value2, function(theta) likelihood(theta))
#plot(Theta_value2, Lik_plot2,xlab = "theta",ylab = "likelihood",
#     main = "Likelihood Function")
#results <- log(Lik_plot2)
#plot(Theta_value2, results)


## Problem 2.b method-of-moments estimator
theta_moment <- asin(mean(x2)-pi)

## Problem 2.c Newton-Raphson with method-of-moments estimator
ObjectFun2 <- function(theta){
  -1*loglik2(theta)  #nlminb is finding minimum value, then multiply by -1
}
results_c <- nlminb(theta_moment,ObjectFun2)

## Problem 2.d different start points -2.7 and 2.7
StartPoints2 <- c(-2.7, 2.7)
results_d <- sapply(StartPoints2, function(x0) nlminb(x0,ObjectFun2))

## Problem 2.e 200 equally spaced starting values between -?? and ??
StartPoints3 <- seq(-pi,pi,length.out = 200)
results_e <- sapply(StartPoints3, function(x0) nlminb(x0,ObjectFun2))
# Partition the values into sets of attraction, That is, divide the set 
# of starting values into separate groups, with each group corresponding 
# to a separate unique outcome of the optimization.

