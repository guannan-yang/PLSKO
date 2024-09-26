#' Simulate Gaussian and binary covariate predictors
#'
#' simulate Gaussian predictors with mean zero
#' and covariance structure determined by "cov_type" argument. Then p_b
#' randomly selected columns are dichotomized.
#'
#' @param n number of observations (rows of X)
#' @param p total number of covariates (columns of X) both continuous and binary
#' @param p_b number of binary covariates (0 <= p_b <= p)
#' @param cov_type character string specifying the covariance function. Can be one of
#' "cov_diag" (independent columns), "cov_equi" (equi-correlated columns), or "cov_ar1" (ar1-correlated columns).
#' The columns are shuffled during simulation
#' @param rho correlation parameter; input to the cov_type function
#'
#' @details This function simulates a data frame, whose rows are multivariate Gaussian with mean zero
#' and covariance structure determined by "cov_type" argument. Then p_b randomly selected columns are
#' dichotomized with the function 1(x>0). The continuous columns are of class "numeric" and the binary
#' columns are set to class "factor".
#'
#' @return the simulated data.frame with n rows and p columns (p_b of which are binary and p-p_b of which are gaussian).
#' Each column is either of class "numeric" or "factor".
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' # all columns are continuous:
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' round(cor(X), 2)
#'
#' # two of the six columns are dichotomized (and set to class factor):
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # The class of each column:
#' unlist(lapply(X, class))
generate_X <- function(n, p, p_b, cov_type, rho=0.5) {

  p_c <- p - p_b

  sigma_z <- do.call(get(cov_type), args=list(n=n, p=p, rho=rho))

  Z <- data.frame(matrix(rnorm(n*p), nrow=n) %*% chol(sigma_z, pivot=TRUE))
  X <- Z

  # Random indices for continuous and binary columns of X
  inds_c <- sample(1:p, size=p_c, replace=FALSE)
  inds_b <- setdiff(1:p, inds_c)

  # Deterministic indices for continuous and binary columns of Z
  seq.c <- 1:p_c
  seq.b <- setdiff(1:p, seq.c)

  # Data frame:
  X[, inds_c] <- identity(Z[, seq.c])
  X[, inds_b] <- lapply((1+sign(Z[, seq.b]))/2, as.factor)

  return(X)

}

# Covariance matrices scaled to be approximately 1/n on the diagonal

cov_diag <- function(n, p, rho=NULL) {
  # Diagonal covariance
  s <- diag(p) / n
  return(s)
}

cov_equi <- function(n, p, rho = 0.5) {
  # Equicorrelated covariance
  s <- (diag(1 - rho, p, p) + rho) / n
  return(s)
}

cov_ar1 <- function(n, p, rho = 0.5) {
  # AR(1) covariance
  s <- toeplitz(rho^(0:(p - 1))) / n
  return(s)
}


#' Simulate Gaussian response from a sparse regression model
#'
#' @param X matrix corresponding to the regression design matrix.
#' Numeric columns of X should have variance = 1/nrow(X), default behavior of generate_X.
#' @param p_nn number of non-null covariate predictors.
#' The regression coefficients (beta) corresponding to
#' columns 1:p_nn of x will be non-zero, all other are set to zero.
#' @param a amplitude of non-null regression coefficients
#'
#' @details This function takes as input data.frame X (created with the function \code{generate_X})
#' that may consist of both numeric and binary factor columns. This data frame is then expanded
#' to a model matrix x (with the model.matrix function). The binary factor variables become dummy
#' indicators that are then scaled by a 0.5*sqrt(nrow(X)) factor so that column-wise variance of the
#' x is equal to 1/n. This makes sense as long as the variance of the numeric
#' columns is also equal to 1/n (which it is if X is generated with the function \code{generate_X}).
#' Next we simulate y ~ N(x%*%beta,I) where the first p_nn beta coefficients are equal to a, while
#' the remaining coefficients (p_nn+1):ncol(x) are set to zero.
#'
#' @return simulated Gaussian response from regression model y = x%*%beta+epsilon, where epsilon ~ N(0,I) and
#' x is the model.matrix of X and the binary dummy indicators of x have been scaled so variance = 1/nrow(X).
#'
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' # Simulate 4 Gaussian and 2 binary covariate predictors:
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # Simulate response from model y = 2*X[,1] + 2*X[,2] + epsilon, where epsilon ~ N(0,1)
#' y <- generate_y(X, p_nn=2, a=2)
generate_y <- function(X, p_nn, a) {

  x <- model.matrix(~., data=X)[,-1]
  which.factor <- as.numeric(which(sapply(X, is.factor)))
  x[,which.factor] <- x[,which.factor]/(0.5*sqrt(nrow(x)))

  n <- nrow(x)
  p <- ncol(x)

  beta <- rep(c(a, 0), c(p_nn, p - p_nn))

  mu <- as.numeric(x %*% beta)

  y <- mu + rnorm(n)

  return(y)

}
