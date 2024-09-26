#' @title Create PC Knockoffs
#' @description Create non-parametric knockoffs based on principal component regression and residuals permutation.
#'
#' @param X An input original design matrix.
#' @param pc.num The number of pricial components to be used for generating knockoff matrices.
#'
#' @return A principal component knockoff matrix.
#'
#' @export
#' @import stats
#'
#' @examples
#' set.seed(10)
#' X <- matrix(rnorm(100), nrow = 10)
#' tmp <- create.pc.knockoff(X = X, pc.num = 5)
create.pc.knockoff <- function(X, pc.num) {
  p <- ncol(X)
  n <- nrow(X)
  bound <- min(n,p)

  if (pc.num >= bound) {
    pc.num <- (bound - 1)
    warning(paste0("Reset num.comp as ", pc.num, "\n"))
  }

  Z <- matrix(NA, ncol = p, nrow = n)
  for (i in 1:p) {
    pca <- stats::prcomp(X[,-i], center = TRUE, scale. = TRUE)
    PCs <- pca$x[, 1:pc.num]

    fit <- stats::.lm.fit(x = PCs, y = X[,i])
    res <- fit$residuals

    tmp <- sample(res) + PCs%*%fit$coefficients
    Z[,i] <- tmp
    X <- cbind(X, tmp)
  }

  return (Z)
}
