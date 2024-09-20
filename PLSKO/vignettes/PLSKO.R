## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#library(PLSKO)

## ----installation-------------------------------------------------------------
# install.packages("devtools")
devtools::install_github("guannan-yang/PLSKO/PLSKO")

library(PLSKO)

## ----eg1_generate-------------------------------------------------------------
set.seed(1234)
n = 100 # number of samples
p = 200 # number of variables ( p > n, high-dimensional setting )
k = 25 # number of important variables with nonzeros coefficients
a = 5 # coefficient for important variables

# generate the variables from a multivariate normal distribution
mu = rep(0, p)
rho = 0.5
Sigma = toeplitz(rho^(0:(p-1))) # we define a covariance matrix with an AR(1) structure
X = MASS::mvrnorm(n, mu, Sigma)

# generate a continuous response variable from a linear model
nonzero = sample(1:p, k)
beta = a * (1:p %in% nonzero) / sqrt(n)
y = X %*% beta + rnorm(n)

# generate a binary response variable from a logistic model
y_bin = rbinom(n, 1, plogis(X %*% beta))

## ----eg1_plsko_pipeline_continous---------------------------------------------
# run the knockoff filter with default settings
result = plsko_filter(X, y) 
print(result)

# compare with the true coefficients
which(beta != 0)
# calculate FDP in this run
fdp <- function(result) {
  sum(beta[result$selected] ==0) / max(length(result$selected), 1)
}
fdp(result)


## -----------------------------------------------------------------------------


