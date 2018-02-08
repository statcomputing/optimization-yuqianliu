beetles <- data.frame(
  days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
  beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))
tt <-  c(8, 28, 41, 63, 69, 97, 117, 135, 154)
N <- c(47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
N0 <- 2

## lodaing package
library(ggplot2)

## Problem 3.a Gauss-Newton approach
#mycode
population <- function(days){
  N0 <- 2
  pop <- K*N0/(N0+(K-N0)*exp(-r*days))
  pop
}
#nls from R
mytab3 <- nls(N ~ K*N0/(N0+(K-N0)*exp(-r*tt)),
              start = list(K = 2000, r = 0.1), trace = TRUE)

## Problem 3.b Countour residual error
#opti: K = 1049.4074581 r = 0.1182683
#range: K from 800 to 1300, r from 0.07 to 0.15
K_plot <- seq(800, 1300, length.out = 100)
r_plot <- seq(0.07, 0.15, length.out = 100)
residual <- rep(0, length(N))
error <- matrix(data = NA, nrow = length(K_plot), ncol = length(r_plot))
residual_error <- function(K, r){
  for (i in 1:length(N)) {
    residual[i] <- (N[i] - (K*N0)/(N0+(K-N0)*exp(-r*tt[i])))^2
  }
  sum(residual)
}
for (i in 1:100) {
  for (j in 1:100) {
    error[j,i] <- residual_error(K_plot[j],r_plot[i])
  }
}

contour(K_plot,r_plot,error)                               
filled.contour(K_plot,r_plot,error,color.palette = terrain.colors,
               xlab="K", ylab="r",
               main="Residual Squared Error")  

## Probelm 3.c