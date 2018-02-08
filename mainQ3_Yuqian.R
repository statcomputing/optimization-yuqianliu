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
              start = list(K = 2000, r = 0.1))
