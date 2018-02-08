## va
x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
StartPoints <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
## lodaing package
library(ggplot2)

## Problem 1.b 
# loglikelihood
Theta_value <- seq(-5, 5, by=0.05)

loglik_cauchy <- function(theta){
  # theta : location
  # x : a vector with samples
  x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
          3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
  n <- length(x1)
  term <- 0
  for(i in 1:n){
    term <- term + log(1+(theta-x1[i])^2)
  }
#  a <- log(1+(theta-x1)^2)
  ll_cauchy <- -n*log(pi)- term
  ll_cauchy
}

LogLik_plot <- sapply(Theta_value, function(theta) loglik_cauchy(theta))
# df <- data.frame(Theta_value, LogLik_plot)
# g <- ggplot(df, aes(Theta_value, LogLik_plot))
# plot log-likelihood
plot(Theta_value, LogLik_plot,xlab = "theta",ylab = "log-likelihood",
     main = "Log-likelihood of Cauchy")

ObjectFun <- function(theta){
  -1*loglik_cauchy(theta)  #nlminb is finding minimum value, then multiply by -1
}

#Newton-Raphson
# z = nlminb(-11, ObjectFun)
mytab <- sapply(c(StartPoints, mean(x1)), function(x0) nlminb(x0,ObjectFun))
dimnames(mytab) <- list(c("theta","objective","convergence","iterations",
                          "evaluations","mesage"),c(StartPoints, mean(x1)))

## Problem 1.c fixed point iterations
#l'(\theta)
logLik_1st <- function(theta){
  x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
          3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
  n <- length(x1)
  for(i in 1:n){
    term <- term + (theta-x1[i])/(1+(theta-x1[i])^2)
  }
  -2*term
}
G <- function(theta,alpha){
  alpha*logLik_1st(theta) + theta
}
myIterator <- function(xinit, f, eps = 1e-06, itmax = 100, verbose = FALSE, 
                       ...) {
  xold <- xinit
  itel <- 1
  repeat {
    xnew <- f(xold, ...)
    if (verbose) 
      cat("Iteration: ", formatC(itel, width = 3, format = "d"), "xold: ", 
          formatC(xold, digits = 8, width = 12, format = "f"), "xnew: ", 
          formatC(xnew, digits = 8, width = 12, format = "f"), "\n")
    if ((supDist(xold, xnew) < eps) || (itel == itmax)) {
      return(xnew)
    }
    xold <- xnew
    itel <- itel + 1
  }
}