## va
x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
StartPoints <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
## lodaing package
library(ggplot2)

## first trial
# Problem 1.b
# loglikelihood
Theta_value <- seq(-5, 5, by=0.01)
term_in_loglik <- function(theta, x){
  # theta : location
  # x: scalar single sample
  log(1+(theta-x)^2)}

loglik_cauchy <- function(theta, X){
  # theta : location
  # x : a vector with samples
  n <- length(X)
  term <- sapply(X, function(x) term_in_loglik(theta, x))
  ll_cauchy <- -n*log(pi)- sum(term)
  ll_cauchy
}

LogLik_plot <- sapply(Theta_value, function(x) loglik_cauchy(x, x1))
df <- data.frame(Theta_value, LogLik_plot)
g <- ggplot(df, aes(Theta_value, l))
 g <- g + geom_line(col='red')
# g <- g + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) #?????????
# g <- g + ggtitle("sigmoid")
g