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

#If warnings of installing dependencies appears, please install the dependencies manually by running the following code:
# install.packages(c("knockoff","progress", "parallel", "doParallel", "foreach"))
#
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("mixOmics")


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

# generate a binary response variable from binomial with sigmoid function as the probability
y_bin = as.factor(rbinom(n, 1, plogis(X %*% beta)))

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


## ----costumised_neighbour-----------------------------------------------------
# define the neighbourhood information based on the AR(1) structure of the variables
# Option 1: define the neighbourhood as a list of length p
nb.list = lapply(1:p, function(i){
  c((i-3):(i-1), (i+1):(i+3))
})
# remove the indices that are out of the range
nb.list = lapply(nb.list, function(x) x[x > 0 & x <= p])

# Then, we run the PLSKO method with the customized neighbourhood information. 
result = plsko_filter(X, y, nb.list = nb.list)
print(result)
fdp(result)

# Option 2: define the neighbourhood as an adjacency matrix
nb.mat = matrix(0, p, p)
for(i in 1:p){
  # make sure the indices are within the range
  nb = (i-3):(i+3)
  nb = nb[nb > 0 & nb <= p]
  nb.mat[i, nb] = 1
}
isSymmetric(nb.mat) # check if the matrix is symmetric

result = plsko_filter(X, y, nb.list = nb.mat)
print(result)
fdp(result)

## ----costumised_ncomp---------------------------------------------------------
# run the PLSKO method with the number of components set to 3 and the sparsity level set to 0.9, which means 90% of the coefficients in PLS regression are zero on each component.
result = plsko_filter(X, y, ncomp = 3, sparsity = 0.9)
print(result)
fdp(result)

## ----binary_response----------------------------------------------------------
# run the knockoff filter with default settings for binary response
result = plsko_filter(X, y_bin)
print(result)
fdp(result)

## ----eg1_plsako_pipeline------------------------------------------------------
# run the PLS-AKO method with default settings
result = plsAKO(X, y)
print(result)
fdp(result)

# Binary response
result = plsAKO(X, y_bin)
print(result)
fdp(result)

## -----------------------------------------------------------------------------


