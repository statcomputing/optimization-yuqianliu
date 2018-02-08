---
title: "Homework2"
author: "Yuqian Liu"
date: "08 February 2018"
documentclass: article
papersize: letter
fontsize: 11pt
output: html_document
---




## Problem 1
**Answer:** 


Cauchy$(\theta ,1)$ distribution has probability density 
$p(x;\theta ) = \frac{1}{{\pi [1 + {{(x - \theta )}^2}]}},$ where
$x \in R$ and $\theta \in R$.

a) The likelihood function of $\theta$ is
$L(x,\theta ) = \prod\limits_{i = 1}^n {p({x_i};\theta )}$.
Then, the log-likelihood is
\begin{align}
\begin{array}{l}
l(\theta ) = \sum\limits_{i = n}^n {\log (p({x_i};\theta ))} \\
 = \sum\limits_{i = n}^n 
 {\log (\frac{1}{{\pi [1 + {{(x - \theta )}^2}]}})} \\
 = \sum\limits_{i = n}^n {[\log \frac{1}{\pi } + 
 \log \frac{1}{{1 + {{(x - \theta )}^2}}}} ]\\
 =  - n\ln \pi  - \sum\limits_{i 
 = n}^n {\ln [1 + {{(\theta  - {x_i})}^2}]} 
\end{array}
\end{align}
Thus, $l'(\theta)$ and $l''(\theta)$ are,
\begin{align}
\begin{array}{l}
l'(\theta ) =  - \sum\limits_{i = n}^n {\frac{{p'({x_i};\theta )}}{{p({x_i};\theta )}}} \\
 =  - \sum\limits_{i = n}^n {[\frac{1}{{1 + 
 {{(\theta  - {x_i})}^2}}} \times 2(\theta  - {x_i})]} \\
 =  - 2\sum\limits_{i = n}^n {\frac{{\theta  - {x_i}}}
 {{1 + {{(\theta  - {x_i})}^2}}}} 
\end{array}
\end{align}

\begin{align}
\begin{array}{l}
l''(\theta ) = \sum\limits_{i = n}^n {\frac{p''(x_i;\theta )}{p(x_i;\theta )}}  
- \sum\limits_{i = n}^n {\{ \frac{p'(x_i;\theta )}{p(x_i;\theta )}} {\} ^2}\\
=  - 2\sum\limits_{i = n}^n {\frac{1 - (\theta  - x_i)^2}{[1 + (\theta  - x_i)^2]^2}} 
\end{array}
\end{align}
The fisher information is
\begin{align}
\begin{array}{l}
I(\theta ) =  - {{\rm E}_\theta }\{ l''(\theta )\} \\
 =  - n{{\rm E}_\theta }[\frac{p''(x_i;\theta )}{p(x_i;\theta )} - 
 {\{ \frac{p'(x_i;\theta )}{p(x_i;\theta )}\} ^2}]\\
 =  - n\int {[\frac{p''(x_i;\theta )}{p(x_i;\theta )} - 
 {{\{ \frac{p'(x_i;\theta )}{p(x_i;\theta )}\} }^2}]p(x_i;\theta )dx} \\
 =  - n\int {\frac{{\{ p'(x)\} }^2}{p(x)}dx} \\
 = \frac{4n}{\pi }\int_{ - \infty }^\infty  {\frac{{{x^2}dx}}{{(1 + x^2)}^3}} \\
 = \frac{n}{2}
\end{array}
\end{align}

b) The observed sample is

```r
x1 <- c(1.77, -0.23, 2.76, 3.80, 3.47, 56.75, -1.34, 4.24, -2.44,
       3.29, 3.71, -2.40, 4.53, -0.07, -1.05, -13.87, -2.53, -1.75)
```
The log-likelihood function is shown in following figure.


```r
plot(Theta_value, LogLik_plot,xlab = "theta",ylab = "log-likelihood",
     main = "Log-likelihood of Cauchy.")
```

<img src="figure/Q1b-1.png" title="plot of chunk Q1b" alt="plot of chunk Q1b" width="90%" style="display: block; margin: auto;" />
Then, Newton-Raphson method with 

```r
StartPoints <- c(-11, -1, 0, 1.5, 4, 4.7, 7, 8, 38)
```
and sample mean as start points have been applied. 
The results are shown in follow.


```r
knitr::kable(mytab1.b, caption = 'Newton-Raphson Method Results', booktabs = TRUE)
```



|    |theta1.b   |itel1.b |
|:---|:----------|:-------|
|par |-0.5914762 |8       |
|par |-0.5914735 |4       |
|par |-0.591474  |5       |
|par |3.021345   |8       |
|par |3.021345   |5       |
|par |3.021346   |7       |
|par |3.021346   |4       |
|par |-0.5914708 |7       |
|par |-0.5914735 |11      |
|par |3.021343   |5       |

c) The fixed-point iterations using $G(x)=\alpha l'(\theta)+\theta$
with scaling choices of $\alpha \in \{1, 0.64, 0.25\}$ is applied
as follow,

```r
# l'(\theta)
logLik_1st <- function(theta, x = x1){
  n <- length(x)
  term <- 0
  for(i in 1:n){
    term <- term + (theta-x1[i])/(1+(theta-x[i])^2)
  }
  -2*term
}
# G = alpha*l'(\theta)+theta
G <- function(theta, alpha, x = x1){
  alpha*logLik_1st(theta, x = x1) + theta
}

myIterator <- function(theta_init, f, eps = 1e-06, itmax = 10000, verbose = FALSE, 
                       ...) {
  theta_old <- theta_init
  itel <- 1
  repeat {
    # theta_t+1 = G
    theta_new <- f(theta_old, ...)
    if (verbose) 
      cat("Iteration: ", formatC(itel, width = 3, format = "d"), "theta_old: ", 
          formatC(theta_old, digits = 8, width = 12, format = "f"), "theta_new: ", 
          formatC(theta_new, digits = 8, width = 12, format = "f"), "\n")
    if ((abs(theta_old - theta_new) < eps) || (itel == itmax)) {
      return(c(theta_new, itel))
    }
    theta_old <- theta_new
    itel <- itel + 1
  }
}
```

The results are shown as follow,.1, .2 and .3 represents results with 3 different 
$\alpha$ values. 

```r
alphas <- c(1, 0.64, 0.25)
```

|    theta.1|    theta.2|    theta.3| iteration.1| iteration.2| iteration.3|
|----------:|----------:|----------:|-----------:|-----------:|-----------:|
| -0.5914740| -0.5914734| -0.5914745|         218|          17|          28|
| -1.1713919| -0.5914734| -0.5914741|        1000|          11|          20|
|  0.1035079| -0.5914736| -0.5914726|        1000|          11|          20|
| -1.1713919|  2.5915177|  3.0213454|        1000|        1000|          12|
|  0.1035079| -0.5914736|  3.0213454|        1000|          12|           8|
|  0.2417269| -0.5914736|  3.0213454|        1000|          13|           9|
|  0.2417269|  3.2398379|  3.0213454|        1000|        1000|           8|
| -1.1063091| -0.5914736|  3.0213454|        1000|          14|          10|
| -1.1063091|  3.2398379|  3.0213454|        1000|        1000|          93|

d) The fisher scoring refined by Newton method has been conducted.
The results are shown as follow.

|      theta| iteration|
|----------:|---------:|
| -0.5914768|        71|
| -0.5914768|        48|
| -0.5914705|        50|
|  3.0213444|        35|
|  3.0213464|        23|
|  3.0213463|        24|
|  3.0213468|        26|
|  3.0213465|        28|
|  3.0213468|       216|


```r
knitr::kable(mytab1.d.withNewton, caption = 'Refined with Newton', booktabs = TRUE)
```

```
## Warning in kable_markdown(x = structure(c("par", "objective",
## "convergence", : The table should have a header (column names)
```



|            |                         |                         |                         |                         |                         |                         |                         |                         |                         |
|:-----------|:------------------------|:------------------------|:------------------------|:------------------------|:------------------------|:------------------------|:------------------------|:------------------------|:------------------------|
|par         |-0.5914768               |-0.5914768               |-0.5914762               |3.021344                 |3.021346                 |3.021346                 |3.021347                 |3.021347                 |3.021347                 |
|objective   |62.73186                 |62.73186                 |62.73186                 |62.59032                 |62.59032                 |62.59032                 |62.59032                 |62.59032                 |62.59032                 |
|convergence |0                        |0                        |0                        |0                        |0                        |0                        |0                        |0                        |0                        |
|iterations  |1                        |1                        |1                        |1                        |1                        |1                        |1                        |1                        |1                        |
|evaluations |2, 1                     |2, 1                     |2, 3                     |2, 1                     |2, 1                     |2, 1                     |2, 1                     |2, 1                     |2, 1                     |
|message     |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |relative convergence (4) |
e) Compare the results from these different methods, we could see
Newton-Raphson method takes general less iterations compared with other methods.
For fixed point method, the choice of $\alpha$ is crucial. We see it may not converge to
the optimal when $\alpha = 1$. However, if we choose a good $\alpha$, the performance are
good.
Fisher Scoring has stable perfomance but with time consuming. Then refined by Newton method,
it get results very fast.

## Problem 2
**Answer:** 


The probability density with parameter $\theta$ is 
$p(x;\theta ){\rm{ }} = {\rm{ }}\frac{{1 - cos(x - \theta )}}{{2\pi }}$,
where $0 \le x \le 2\pi$ and $\theta  \in ( - \pi ,{\rm{ }}\pi )$.

a) The log-likelihood function of $\theta$ is
\begin{align}
\begin{array}{l}
l(\theta ) = \sum\limits_{i = n}^n {\log (p({x_i};\theta ))} \\
 = - 2n\ln \pi  + \sum\limits_{i = n}^n {\ln [1 - cos(x - \theta )]} 
\end{array}
\end{align}
The observed sample is

```r
x2 <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
```
The log-likelihood function is shown in Figure.

<img src="figure/Q2a-1.png" title="plot of chunk Q2a" alt="plot of chunk Q2a" width="90%" style="display: block; margin: auto;" />

b) Method of moments
\begin{align}
\begin{array}{l}
{\rm E}(X|\theta ) = \int_0^{2\pi } {p(x;\theta )} xdx\\
 =  -\frac{1}{2\pi }[\cos (x - \theta ) - \frac{x^2}{2} + x\sin (x - \theta )]|_0^{2\pi }\\
 = \pi  + \sin \theta 
\end{array}
\end{align}
Thus, ${\hat \theta _{moment}}$ could be calculated based on ${\rm E}(X|\theta ) = \bar x$. 
Solve for $\theta$,
\begin{align}
${\hat \theta _{moment}} = \arcsin (\bar x - \pi )$
\end{align}

c) With $\theta_0 = \hat \theta _{moment}$, using the Newton-Raphson method the MLE is

```
## Iteration:    4 ,theta:    0.00311816 ,objective:   38.67500474
```

d) With $\theta_0 = -2.7$ and $\theta_0 = 2.7$, the results are

```
## Warning in kable_markdown(x = structure(c("par", "objective",
## "convergence", : The table should have a header (column names)
```



|            |                         |                         |
|:-----------|:------------------------|:------------------------|
|par         |-1.662712                |2.848415                 |
|objective   |56.68446                 |67.80782                 |
|convergence |0                        |0                        |
|iterations  |6                        |5                        |
|evaluations |8, 9                     |7, 7                     |
|message     |relative convergence (4) |relative convergence (4) |

e) Repeat the above using 200 equally spaced starting values between $-\pi$ and $\pi$.
The results are too large to present here. It is seperated by sets of attraction.

```r
StartPoints3 <- seq(-pi,pi,length.out = 200)
results_e <- sapply(StartPoints3, function(x0) nlminb(x0,ObjectFun2))
```

## Problem 3
**Answer:** 


The model for population growth is given by $\frac{{dN}}{{dt}} = rN(1 - \frac{N}{K})$.
The solution of this is 
${N_t} = f(t){\rm{ }} = \frac{{K{N_0}}}{{{N_0} + (K - {N_0}){\rm{ }}exp( - rt)}}$

a) By using nls() in R, the results are K = 1049.4074581, r = 0.1182683 at optimal.


Code chunk:

```r
beetles <- data.frame(
  days = c(0, 8, 28, 41, 63, 69, 97, 117, 135, 154),
  beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))
tt <-  c(8, 28, 41, 63, 69, 97, 117, 135, 154)
N <- c(47, 192, 256, 768, 896, 1120, 896, 1184, 1024)
N0 <- 2

## lodaing package
library(ggplot2)

## Problem 3.a Gauss-Newton approach
#nls from R
mytab3 <- nls(N ~ K*N0/(N0+(K-N0)*exp(-r*tt)),
              start = list(K = 2000, r = 0.1), trace = TRUE)
```

b) The contour plots of the sum of squared errors are shown as follow,
<img src="figure/Q3b1-1.png" title="plot of chunk Q3b1" alt="plot of chunk Q3b1" width="90%" style="display: block; margin: auto;" />
<img src="figure/Q3b2-1.png" title="plot of chunk Q3b2" alt="plot of chunk Q3b2" width="90%" style="display: block; margin: auto;" />

Code chunk:

```r
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

contour(K_plot,r_plot,error, xlab="K", ylab="r",
               main="Residual Squared Error")                               
filled.contour(K_plot,r_plot,error,color.palette = terrain.colors,
               xlab="K", ylab="r",
               main="Residual Squared Error")  
```
