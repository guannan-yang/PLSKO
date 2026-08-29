# This is a function from "knockoff"
# that simply standardize the design matrix
# so that column sum of squares = 1 (instead of n - 1 from "scale")
normc <- function(X, center = T)
{
  X.centered = scale(X, center = center, scale = F)
  X.scaled = scale(X.centered, center = F, scale = sqrt(colSums(X.centered^2)))
  X.scaled[, ]
}

#' Generating multiple knockoff statistics using coefficient estimate from glmnet
#'
#' This function computes the knockoff statistics based on the absolute value of the coefficient estimate from glmnet.
#'
#' @param X A \code{n}-by-\code{p} matrix of original variables
#' @param X_k A \code{n}-by-\code{d} matrix of multiple knockoff variables
#' @param y A \code{n} vector of response
#' @param omega A \code{p} vector indicating the weights for each variable. For now, we require each entry to be an integer greater than or equal to 2.
#' @param family The conditional distribution of y given X. See the family option for \code{glmnet}.
#' @param nlam Number of tuning parameter lambda used in fitting the lasso. Default to be 500.
#' @param lam_min_ratio The ratio of the minimum and the maximum value of lambda in constructing the tuning parameters. Default to be \code{1e-4}.
#'
#' @import glmnet
#'
#' @return An list of three components:
#' \describe{
#' \item{\code{kappa}}{the vector of indices of winner for each variable competing with its multiple knockoff counterparts. \code{kappa[j] = 1} indicates that the original variable is beating all of its knockoff counterparts, and \code{kappa[j]} not equal to 1 means otherwise.}
#' \item{\code{tau}}{a vector of scores determining the order for which we consider to include variables into the model.}
#' \item{\code{score_total}}{the matrix containing the original `glmnet` coefficient estimates for each variable and its knockoff counterparts. For example, \code{score_total[1:omega[j], j]} is the coefficients estimates for the j-th variables and its \code{omega_j} - 1 knockoff counterparts.}
#' }
#'
#' @examples
#' library(cheapknockoff)
#' set.seed(123)
#' n <- 100
#' p <- 30
#' x <- matrix(data = rnorm(n * p), nrow = n, ncol = p)
#' y <- x[, 1] - 2 * x[, 2] + rnorm(n)
#' omega <- c(2, 9, sample(seq(2, 9), size = 28, replace = TRUE))
#' # construct multiple knockoffs
#' X_k <- multiple_knockoff_Gaussian(X = x, mu = rep(0, p), Sigma = diag(1, p), omega = omega)
#' # compute knockoff statistics
#' stat <- cheapknockoff::stat_glmnet_coef(X = x, X_k = X_k, y = y, omega = omega)
#' @export
stat_glmnet_coef <- function(X, X_k, y, omega, family = "gaussian",
                             nlam = 500, lam_min_ratio = 5e-4){
  # this function adapts largely from "stat.lasso_coefdiff" function
  # from the package "knockoff"

  # return the absolute value of the glmnet estimate of coefficient
  # as the importance score
  p <- ncol(X)

  # here we assume that for each variable we construct same number of knockoffs
  # first combine original variables with the knockoffs
  # x <- scale(cbind(X, X_k))
  x <- normc(cbind(X, X_k))
  # problem dimensions
  n <- nrow(x)
  d <- ncol(x)

  # fit a glmnet with cross-validation
  fit <- glmnet::cv.glmnet(x, y, family = family, intercept = TRUE,
                           nlambda = nlam,
                           lambda.min.ratio = lam_min_ratio,
                           standardize = FALSE, nfolds = 5)
  # here the importance score is the estimated coefficients
  # excluding the intercept
  Z <- coef(fit, s = "lambda.min")[2:(d + 1)]

  # Tau is a (kappa + 1)-by-p matrix
  # with T[i, j] the importance score of the i-th knockoff of the j-th variable
  score_total <- matrix(abs(Z), ncol = p, byrow = TRUE)

  kappa <- rep(NA, p)
  tau <- rep(NA, p)
  for(i in seq(p)){
    score <- score_total[1:omega[i], i] # when all omega equals to (n_ko+1), this becomes the same as score_total[, i]
    if(all(score == 0)){
      kappa[i] <- 2
      tau[i] <- 0
    }
    else{
      # kappa[i] is the index of the largest score
      kappa[i] <- which.max(score)
      # tau[i] is the difference between the largest and the second largest score
      tau[i] <- (max(score) - max(score[-which.max(score)])) / omega[i]
      # tau[i] <- max(score) - max(score[-which.max(score)])
      # tau[i] <- (max(score) - min(score)) / omega[i]
    }
  }
  return(list(kappa = kappa, tau = tau, score_total = score_total))
}
