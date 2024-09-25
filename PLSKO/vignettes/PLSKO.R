## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#library(PLSKO)

## ----installation-------------------------------------------------------------
# install.packages("devtools")
#devtools::install_github("guannan-yang/PLSKO/PLSKO", quiet = TRUE, upgrade = "never")

library(PLSKO)

#If warnings of installing dependencies appears, please install the dependencies manually by running the following code:
# install.packages(c("knockoff","progress", "parallel", "doParallel", "foreach"))
#
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("mixOmics")


## ----eg1_generate-------------------------------------------------------------
set.seed(1)
n = 100 # number of samples
p = 150 # number of variables ( p > n, high-dimensional setting )
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
  if(class(result) == "knockoff.result") fdp = sum(beta[result$selected] ==0) / max(length(result$selected), 1)
  else if (class(result) == "AKO.result") fdp = sum(beta[result$ako.s] ==0) / max(length(result$ako.s), 1)
  return(fdp)
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
result = plsko_filter(X, y, ncomp = 3, sparsity = 0.95)
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

#Binary response
# result = plsAKO(X, y_bin)
# print(result)
# fdp(result)

## ----eg2_semi_synthetic_generate----------------------------------------------
data("prot_placenta")
X = as.matrix(prot_placenta$abundance)

#generate the response variable y from a linear model
set.seed(1)
n = nrow(X)
p = ncol(X)
k = 8 # number of important variables with nonzeros coefficients
nonzero = sample(1:p, k) # randomly select 8 important variables

beta = as.numeric(1:p %in% nonzero) # assign non-zero coefficients to the important variables
y = X %*% beta 


## ----eg2_plsko_pipeline-------------------------------------------------------
result = plsko_filter(X, y) 
print(result)
fdp(result)

## ----eg2_plsko_pipeline_custom------------------------------------------------
result <- plsko_filter(X, y, threshold.abs = 0, ncomp = 5, sparsity = 0.8) # set the absolute correlation threshold to 0 (every one is neighour with each other), the number of components to 5, and the sparsity level in PLS regression to 0.8
print(result)
fdp(result)

## ----eg2_plsako_pipeline------------------------------------------------------
result = plsAKO(X, y, threshold.abs = 0, ncomp = 5, sparsity = 0.8)
print(result)
fdp(result)

# average the results from 25 iterations
mean(unlist(lapply(result$s, function(x) fdp(x))))

## ----eg3_advanced_usage_generatedata------------------------------------------
data("cfRNA_placenta")
X = as.matrix(cfRNA_placenta$counts)
#generate the response variable y from a linear model
set.seed(1)
n = nrow(X)
p = ncol(X)
k = 8 # number of important variables with nonzeros coefficients
nonzero = sample(1:p, k) # randomly select 8 important variables

beta = as.numeric(1:p %in% nonzero) # assign non-zero coefficients to the important variables
y = X %*% beta 

## ----eg3_advanced_usage_knockoff----------------------------------------------
# generate the knockoff variables with default settings of the plsko function
plsko_default = plsko(X) 

# generate the knockoff variables with customized settings
plsko_custom = plsko(X, threshold.abs = 0, ncomp = 7, sparsity = 0.8)

# generate the knockoff variables with customized settings and neighbourhood information
# Here we may estimate the neibourhood information from the data by using graphical lasso, which provides a sparse estimate of the precision matrix
library(glasso)
cov.mat = cov(X)
glasso_res = glasso(cov.mat, rho = 0.1) # estimate the precision matrix using graphical lasso
nb.list = lapply(1:p, function(i) which(glasso_res$wi[i, ] != 0)) # get the neighbourhood information from the precision matrix
plsko_custom_nb = plsko(X, ncomp = 7, nb.list = nb.list)


# Or you can use other knockoff methods to generate the knockoff variables
# For example, you can use the method provided in the `knockoff` package to generate the knockoff variables
library(knockoff)
ko_soa = knockoff::create.second_order(X)

## ----eg3_advanced_usage_ko_filter---------------------------------------------
# calculate the importance scores and perform the knockoff filtering and variable selection with the default settings
result_default = ko_filter(X, plsko_default, y)
print(result_default)
fdp(result_default)

# calculate the importance scores and perform the knockoff filtering and variable selection with the customized settings
result_custom = ko_filter(X, plsko_custom, y)
print(result_custom)
fdp(result_custom)

# calculate the importance scores and perform the knockoff filtering and variable selection with the customized settings and neighbourhood information
result_custom_nb = ko_filter(X, plsko_custom_nb, y)
print(result_custom_nb)
fdp(result_custom_nb)

# calculate the importance scores and perform the knockoff filtering and variable selection with the knockoff variables generated by the `knockoff` package
result_ko_soa = ko_filter(X, ko_soa, y)
print(result_ko_soa)
fdp(result_ko_soa)

# Or you might have different options to calculate the importance scores
result_default_rf <- ko_filter(X, plsko_default, y, method = "RF") # e.g. use random forest to calculate the importance scores
print(result_default_rf)
fdp(result_default_rf)


## ----eg3_advanced_usage_ko_withW----------------------------------------------
# calculate the importance scores using self-defined method, e.g. difference of absolute value marginal correlation between X and y between the original and knockoff variables 

my_knockoff_stat = function(X, X_k, y) {
  abs(t(X) %*% y) - abs(t(X_k) %*% y)
}
W = my_knockoff_stat(X, plsko_custom, y)
result_my = ko_withW(W, q = 0.05)
print(result_my)

## ----eg4_plsako_steps_1-------------------------------------------------------
X = as.matrix(cfRNA_placenta$counts)
y = cfRNA_placenta$metadata$PE # the real binary response variable indicates pre-eclampsia or control

# generate multiple knockoff variables by PLSKO with customised setting 
n_ko = 15
plsko_list = lapply(1:n_ko, function(i) plsko(X, seed = i)) # generate 15 knockoff independently. Note that seed needs to be set since the default seed is 1 in the function, without specifying, the same knockoff variables will be generated in each iteration.


## ----eg4_plsako_steps_2-------------------------------------------------------
# calculate the importance scores and perform the knockoff filtering and variable selection with the PLS-AKO method
result_ako = AKO_withKO(X, plsko_list, y)
print(result_ako)

## ----eg4_plsako_steps_selfKO--------------------------------------------------
soa_list = lapply(1:n_ko, function(i) knockoff::create.second_order(X))
result_ako_soa = AKO_withKO(X, soa_list, y)
print(result_ako_soa)

## ----eg4_plsako_steps_3-------------------------------------------------------
# calculate the importance scores from multiple knockoff variables using self-defined method
# We use the difference of absolute value of Z-score (on contrast of PE or control group) between the original and knockoff variables as the importance score
my_knockoff_stat_Z <- function(X, Xk, y){
  X_new <- cbind(X, Xk)
 beta <- coef(lm(X_new~y))[2,] # run orgianal variables and knockoff variables together into an OLS regression and extract Z-score
 abs(beta[1:ncol(X)]) - abs(beta[(ncol(X)+1):ncol(X_new)])
}

# calculate the importance scores from multiple knockoff variables
W_list = lapply(1:n_ko, function(i) my_knockoff_stat_Z(X, plsko_list[[i]], y))

# perform the knockoff filtering and variable selection with AKO method
result_ako_W = AKO_withW(W_list, q = 0.05)
print(result_ako_W)

