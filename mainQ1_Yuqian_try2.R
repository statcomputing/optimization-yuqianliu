## va
x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
StartPoints <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
alphas <- c(1, 0.64, 0.25)

## lodaing package
library(ggplot2)

## Problem 1.b 
# loglikelihood
Theta_value <- seq(-5, 5, by=0.05)

loglik_cauchy <- function(theta, x = x1){
  # theta : location
  # x : a vector with samples
  n <- length(x)
  term <- 0
  for(i in 1:n){
    term <- term + log(1+(theta-x[i])^2)
  }
#  a <- log(1+(theta-x1)^2)
  ll_cauchy <- -n*log(pi)- term
  ll_cauchy
}

LogLik_plot <- sapply(Theta_value, function(theta) loglik_cauchy(theta, x = x1))
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
# colnames(mytab) <- c(StartPoints, mean(x1))
# dimnames(mytab) <- list(c("theta","objective","convergence","iterations",
#                          "evaluations","mesage"),c(StartPoints, mean(x1)))
theta1.b <-  mytab[1,1]
itel1.b <-  mytab[4,1]
for (i in 1:9) {
  theta1.b <-  c(theta1.b, mytab[1,i+1])
  itel1.b <-  c(itel1.b, mytab[4,i+1])
}
mytab1.b <- cbind(theta1.b,itel1.b)
# names(theta1.b) <- list(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38, "mean(x1)")
# names(itel1.b) <- list(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38, "mean(x1)")
#mytab1.b <- data.frame(`Starts` = c(StartPoints, mean(x1)) , `Theta`= theta1.b, 
#                       `Iteration` = itel1.b, stringsAsFactors = FALSE, check.names = FALSE)

## Problem 1.c fixed point iterations
#l'(\theta)
logLik_1st <- function(theta, x = x1){
  n <- length(x)
  term <- 0
  for(i in 1:n){
    term <- term + (theta-x1[i])/(1+(theta-x[i])^2)
  }
  -2*term
}

# G = alpha*l'(theta)+theta
G <- function(theta, alpha, x = x1){
  alpha*logLik_1st(theta, x = x1) + theta
}

myIterator <- function(theta_init, f, eps = 1e-06, itmax = 1000, verbose = FALSE, 
                       ...) {
  theta_old <- theta_init
  itel <- 1
  repeat {
    #theta_t+1 = G
    theta_new <- f(theta_old, ...)
    #if (verbose) 
    #  cat("Iteration: ", formatC(itel, width = 3, format = "d"), "theta_old: ", 
    #      formatC(theta_old, digits = 8, width = 12, format = "f"), "theta_new: ", 
    #      formatC(theta_new, digits = 8, width = 12, format = "f"), "\n")
    if ((abs(theta_old - theta_new) < eps) || (itel == itmax)) {
      return(c(theta_new, itel))
    }
    theta_old <- theta_new
    itel <- itel + 1
  }
}
# test
# sy <- myIterator(-11, G, eps = 1e-06, itmax = 10000, verbose = TRUE, alpha = 1)
mytab1.c_theta <- matrix(NA, length(StartPoints), length(alphas))
mytab1.c_itel <- matrix(NA, length(StartPoints), length(alphas))

for (i in 1:length(StartPoints)) {
  theta_init <- StartPoints[i]
  for (j in 1:length(alphas)) {
    alpha <- alphas[j]
     sy <- myIterator(theta_init, G, eps = 1e-06, itmax = 1000, verbose = TRUE, alpha,
                      x = x1)
     mytab1.c_theta[i,j] <- sy[1]
     mytab1.c_itel[i,j] <- sy[2]
  }
}
mytab1.c <- data.frame('theta' = mytab1.c_theta, 'iteration' = mytab1.c_itel,
                       stringsAsFactors = FALSE, check.names = FALSE )
## Problem 1.d Fisher Scoring + Newton
#Fisher scoring
I <- length(x1)/2
mytab1.d_theta <- rep(NA, length(StartPoints))
mytab1.d_itel <- rep(NA, length(StartPoints))
myFisher <- function(theta_init, f, fisher_value, eps = 1e-06, itmax = 1000, verbose = FALSE, 
                     ...) {
  theta_old <- theta_init
  itel <- 1
  repeat {
    #theta_t+1 = theta_t + l'(theta)/I
    theta_new <- theta_old + f(theta_old, ...)/fisher_value
    if ((abs(theta_old - theta_new) < eps) || (itel == itmax)) {
      return(c(theta_new, itel))
    }
    theta_old <- theta_new
    itel <- itel + 1
  }
}

for (i in 1:length(StartPoints)) {
  theta_initd <- StartPoints[i]
  sy <- myFisher(theta_initd, logLik_1st, I, eps = 1e-06, itmax = 1000, verbose = TRUE,
                   x = x1)
  mytab1.d_theta[i] <- sy[1]
  mytab1.d_itel[i] <- sy[2]
}
mytab1.d <- data.frame('theta' = mytab1.d_theta, 'iteration' = mytab1.d_itel,
                       stringsAsFactors = FALSE, check.names = FALSE )
#Newton
mytab1.d.withNewton <- sapply(mytab1.d_theta, function(x0) nlminb(x0,ObjectFun))

