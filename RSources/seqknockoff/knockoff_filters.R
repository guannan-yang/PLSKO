#' The Knockoff Filter
#'
#' This is an adaptation of the knockoff.filter function from the R-package knockoff.
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial").
#' @param fdr target false discovery rate. Can be a vector of multiple thresholds.
#' @param family should be "gaussian" if y is numeric, but "binomial" if y is a binary factor variable.
#' @param knockoffs user-specified function to construct knockoff of X. It must take as input
#' the data.frame (or tibble) X and return a data.frame (or tibble) X_k of corresponding
#' knockoffs. By default, \code{knockoffs=knockoffs_seq}, but other option is
#' \code{knockoffs=knockoffs_mx} for X with numeric columns only (see ?knockoffs_seq and ?knockoffs_mx).
#' @param statistic user-specified function that constructs feature statistics used to assess variable importance.
#' By default \code{statistic=stat_glmnet} (see ?stat_glmnet).
#'
#' @details This function takes input X with either "numeric" or "factor" columns (or both), input y can be either
#' numeric or binary factor, and user may specify multiple fdr thresholds. The function performs the
#' knockoff filter, which consists of three steps: 1) Simulate knockoff of X with the input function \code{knockoffs},
#' 2) calculate importance statistic W with the input function \code{statistic}, and 3) calculate the knockoff+ threshold
#' for each target fdr provided. Finally, selection is made based on W > threshold.
#'
#' @return if length(fdr)=1 the function returns a vector of selected indices, otherwise a list of selected indices, one selection vector per fdr threshold supplied.
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n=1000, p=10, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' # Simulate response from model y = X%*%beta + epsilon, where epsilon ~ N(0,1) with
#' # first 5 beta-coefficients = 8 (all other zero).
#' y <- generate_y(X, p_nn=5, a=8)
#'
#' S <- knockoff_filter(X, y, fdr=c(0.05, 0.1, 0.2))
#'
#' # dichotomize y for logistic regression knockoff filter:
#' y <- factor(y > median(y))
#'
#' # Below the family argument gets passed to the statistic = knockoff::stat.glmnet_coefdiff function:
#' S <- knockoff_filter(X, y, fdr=c(0.05, 0.1, 0.2), family="binomial")
knockoff_filter <- function(X, y, fdr = 0.2, family="gaussian", knockoffs = knockoffs_seq, statistic = stat_glmnet) {

  check_design(X)
  check_family(y, family)

  # Calculate the knockoff copy of X:
  X_k <- knockoffs(X)

  W = statistic(X=X, X_k=X_k, y=y, family=family)

  if (length(fdr)==1) {
    ko.thres <- knockoff::knockoff.threshold(W, fdr=fdr, offset=1)
    selected <- which(W >= ko.thres)
  } else {
    ko.thres <- lapply(fdr, function(x) knockoff::knockoff.threshold(W, fdr=x, offset=1))
    selected <- lapply(ko.thres, function(x) which(W >= x))
  }

  return(selected)

}


#' Importance statistics: Absolute elastic-net coefficient differences between original and knockoff variables
#'
#' This is a wrapper function that calls knockoff::glmnet.stat_coefdiff. The input data.frames
#' are first converted to design matrices (with the function model.matrix). This means that if
#' the input features contain factor variables then an importance statistic is calculated for
#' each dummy variable as determined by the model.matrix contrasts (defaults to indicator dummy
#' variables with a reference level).
#'
#' @param X original data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param X_k knockoff data.frame (or tibble) with "numeric" and "factor" columns only. The dimensions and column classes must match
#' those of the original X.
#' @param y response vector with \code{length(y) = nrow(X)}. Accepts "numeric" (family="gaussian") or binary "factor" (family="binomial").
#' @param family should be "gaussian" if y is numeric, but "binomial" if y is a binary factor variable.
#' @param ... additional parameters passed to knockoff::stat.glmnet_coefdiff
#'
#' @return a vector of importance statistics W of length equal to number of columns of the model.matrix of X.
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' # Simulate 10 Gaussian covariate predictors:
#' X <- generate_X(n=1000, p=10, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#'  # Calculate the knockoff copy of X:
#'  X_k <- knockoffs(X)
#'
#' # Simulate response from model y = X%*%beta + epsilon, where epsilon ~ N(0,1) with
#' # first 5 beta-coefficients = 8 (all other zero).
#' y <- generate_y(X, p_nn=5, a=8)
#'
#' W <- stat_glmnet(X=X, X_k=X_k, y=y, family="gaussian", nfolds=5)
stat_glmnet <- function(X, X_k, y, family, ...) {

  # Calculate the W-statistic from LASSO:
  x <- model.matrix(~., data=X)[,-1]
  x_k <- model.matrix(~., data=X_k)[,-1]

  W <- knockoff::stat.glmnet_coefdiff(X=x, X_k=x_k, y=y, family=family, ...)

  return(W)

}


#' Select variables based on the heuristic multiple selection algorithm from Kormaksson et al. 'Sequential
#' knockoffs for continuous and categorical predictors: With application to a large psoriatic arthritis clinical
#' trial pool.' Statistics in Medicine. 2021;1–16.
#'
#' @param S list of variable selection indices
#' @param p number of variables. Each element of the list of selection indices should be a subset of 1:p.
#' @param trim trimming probability threshold. A sensible default is \code{trim=0.5}.
#'
#' @return a single "most frequent" variable selection among the multiple selections in S.
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' S <- list(c(sample(1:10,7), sample(11:20,1)),
#' c(sample(1:10,7), sample(11:20,1)),
#' c(sample(1:10,7), sample(11:20,1)))
#'
#' multi_select(S, p=20)
multi_select <- function(S, p, trim=0.5) {

  countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }

  candidates <- 1:p

  Nknockoff <- length(S)

  var.freq.table <- table(factor(unlist(S), levels=candidates))

  trim.seq <- unique(var.freq.table/Nknockoff)
  trim.seq <- trim.seq[trim.seq >= trim]

  if (length(trim.seq) > 0) {
    selected <- lapply(trim.seq, function(x) find_single_optimal_variable_set(S, p, trim=x))
    selected <- selected[[which.max(unlist(lapply(selected, length)))]]
  } else {
    selected <- integer(0)
  }

  return(selected)

}


