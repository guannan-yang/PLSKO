#' Sequential knockoffs for continuous and categorical variables
#'
#' @param X data.frame (or tibble) with "numeric" and "factor" columns only. The number of columns, ncol(X) needs to be > 2.
#' @param seq_simulator function that simulates sequential knockoffs. Default is the function \code{sim_EN}, which simulates response from an estimated elastic-net model
#' @param ... other parameters passed to the function seq_simulator. For the default (elastic-net sequential seq_simulator, \code{seq_simulator = sim_EN})
#' these other parameters are passed to cv.glmnet.
#'
#' @details \code{knockoffs_seq} performs sequential knockoff simulation using elastic-net regression.
#' @return sequential knockoff copy of X. A data.frame or tibble of same type and dimensions as X.
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=2, cov_type="cov_equi", rho=0.5)
#'
#' # knockoffs based on sequential elastic-net regression with penalty alpha:
#' Xk <- knockoffs_seq(X, alpha=0.5)
knockoffs_seq <- function(X, seq_simulator = sim_EN, ...) {

  check_design(X)

  knockoffs <- X

  # add ".tilde" to column names:
  names(knockoffs) <- paste0(names(knockoffs),".tilde")

  # Randomly shuffle column indices of X:
  shf <- sample(ncol(X))

  # Loop through the columns of input data (in random order)
  loop.count <- 1
  for (i in shf) {

    y <- X[[i]] # i-th column serves as response
    Xp <- X[,-i] # columns[-i] serve as predictors

    if (loop.count > 1) Xp <- cbind(knockoffs[,shf[1:(loop.count-1)]], Xp)

    knockoffs[[i]] <- seq_simulator(y = y, X = Xp, ...)

    loop.count <- loop.count + 1

  }

  # remove ".tilde" from column names:
  names(knockoffs) <- gsub(".tilde","", names(knockoffs))

  return(knockoffs)

}


#' Simulate from elastic-net regression model
#'
#' @param y response vector (either "numeric" or "factor") that gets passed to cv.glmnet
#' @param X data.frame of covariates that are passed to cv.glmnet
#' @param ... other parameters passed to the function cv.glmnet
#'
#' @return simulated response
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' X = data.frame(matrix(rnorm(100 * 20), 100, 20))
#' y = X[,1] + rnorm(100)
#'
#' # simulate from elastic-net regression with elastic-net penalty alpha:
#' ysim = sim_EN(y=y, X=X, alpha=0.5)
#'
#' # simulated versus input response:
#' plot(y, ysim)
sim_EN <- function(y, X, ...) {

  x <- model.matrix(~., data = X)[,-1]

  if (is.factor(y)) {

    classes <- levels(y)

    K <- length(classes)

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="multinomial", intercept=TRUE, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min")[[2]])[-1]

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")

    mat.multinom <- apply(mu, 1, function(prob) rmultinom(n=1, size=1, prob=prob))

    y.sim <- classes[apply((1:K)*mat.multinom, 2, max)]

    y.sim <- factor(y.sim, levels=classes)

    rmse <- NULL

  } else {

    if(!is.numeric(y)) stop("class(y) needs to be either 'numeric' or 'factor'")

    gm.cv <- glmnet::cv.glmnet(y=y, x=x, family="gaussian", intercept=TRUE, ...)

    # Beta coefficients (excluding intercept)
    beta.coefs <- as.numeric(coef(gm.cv, s = "lambda.min"))[-1]

    # columns of predictor matrix corresponding to non-zero beta.coefs:
    non.zero.cols <- which(beta.coefs != 0)

    # Total number of non-zero parameters (including intercept, hence + 1)
    s.lambda = length(non.zero.cols) + 1

    mu <- predict(gm.cv, newx=x, type="response", s="lambda.min")
    
    #rmse = sqrt(sum((y-mu)^2)/(length(y) - s.lambda))  # NANs might be produce when n << p
    rmse = sqrt(sum((y-mu)^2)/max((length(y) - s.lambda), 1))

    y.sim <- rnorm(n=length(y), mean=mu, sd=rmse)

  }

  return(y.sim)

}


#' Gaussian MX-knockoffs for continuous variables
#'
#' @param X  data.frame (or tibble) with "numeric" columns only. The number of columns, ncol(X) needs to be > 2.
#'
#' @details \code{knockoffs_mx} performs MX knockoff simulation.
#'
#' @return Second-order multivariate Gaussian knockoff copy of X
#' @export
#'
#' @examples
#' #' library(seqknockoff)
#'
#' set.seed(1)
#'
#' X <- generate_X(n=100, p=6, p_b=0, cov_type="cov_equi", rho=0.5)
#'
#' Xk <- knockoffs_mx(X)
knockoffs_mx <- function(X) {

  check_design(X, method="mx")

  knockoffs <- data.frame(knockoff::create.second_order(as.matrix(X)))

  return(knockoffs)

}
