#' Multiple knockoff filter by Gimenez & Zou (2018)
#'
#' This function implements the stable m-knockoff filter by Gimenez & Zou (2018)
#' given the test statistics (kappa, tau)
#'
#' @param kappa A \code{p} vector of test statistics, with kappa_i = 1 indicating the original variable winning
#' @param tau A \code{p} vector of test statistics, showing the manitude/importance of the variable
#' @param fdr The pre-specified upper bound level of controlled FDR.
#' @param n_knockoff number of knockoffs constructed
#'
#' @return The selected variable and their corresponding scores
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
#' stat <- stat_glmnet_coef(X = x, X_k = X_k, y = y, omega = omega)
#' # run gz filter
#' gz_output <- filter_gz(kappa = stat$kappa, tau = stat$tau, fdr = 0.2, n_knockoff = max(omega))
#' @export

filter_gz <- function(kappa, tau, fdr = 0.2, n_knockoff, offset = 1){
  m <- length(kappa)
  #stopifnot(n_knockoff > 1)

  # input check
  stopifnot(length(kappa) == length(tau))
  stopifnot(all(tau >= 0))

  # now tau[inc] is in non-increasing order
  inc <- order(-tau)

  # kap is the permutations of kappa, in the order inc
  kap <- kappa[inc]

  # now we start computing the threshold:
  # assume that tau is non-increasing
  # if there is no discovery at all...
  # return k = 0
  if(all(kap > 1))
    score <- 0
  # the threshold is the ratio of the following two things
  # num[i] = (1 + \sum_j \indi {j \leq i: W_j < 0}) / num_knockoff
  num <- (cumsum((kap > 1)) + offset) / n_knockoff
  # den[i] = \sum_j \indi {j \leq i: W_j > 0}
  den <- cumsum((kappa == 1))

  score <- num / den
  score[den == 0] <- 0

  # now we compute \hat{k} and finish the filter
  if(all(score > fdr))
    hatk <- 0
  else
    hatk <- max(which(score <= fdr))

  # the discovery set S in the order inc
  hatS <- which(kap == 1)
  hatS <- hatS[hatS <= hatk]
  # now we need to go back to the original order
  selected <- inc[hatS]
  return(list(selected = sort(selected), score = score, tau = tau))
}
