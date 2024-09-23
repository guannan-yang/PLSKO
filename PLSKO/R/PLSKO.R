#' Generate a knockoff variable set with PLSKO using PLS regression
#'
#' @importFrom mixOmics spls
#' @import progress
#' @param X A numeric matrix or data frame. The original design data matrix with \eqn{n} observations as rows and \eqn{p} variables as columns.
#' @param nb.list Optional. A list of length \eqn{p} or adjacency matrix of \eqn{p \times p} that defines the neighbourship of variables. A list of length \eqn{p} should include the neighbours' index of each variable from \eqn{X_1} to \eqn{X_p} in order; The \eqn{i^{th}} element in the list includes the indices of the neighbour variables of \eqn{X_i}, or \code{NULL} when no neighbours. A adjacency matrix should be symmetric with only binary element and  where \eqn{M_{ij} = 1} when \eqn{X_i} and \eqn{X_j} is defined as neighbours; otherwise \eqn{M_{ij} = 0} when not neighbour or on diagonal (i.e. \eqn{i = j}). If not provided or NULL, the neighborhoods are determined based on correlations.
#' @param threshold.abs Optional. Used when \code{nb.list} is not provided. A value between \eqn{0} and \eqn{1}. A numeric value specifying an absolute value of correlation threshold to define neighborhoods. If not provided, the function uses \code{threshold.q}.
#' @param threshold.q Optional. Used when \code{nb.list} and \code{threshold.q} are both not provided. A numeric value between 0 and 1 indicating the quantile of the sample correlation values to use as a threshold. When not provided, the function uses the 0.80 quantile.
#' @param ncomp Optional. An integer specifying the number of components to use in the PLS regression. Default is NULL, the \code{ncomp} is determined empirically by \eqn{PC_p1} criterion.
#' @param sparsity Optional. A numeric value between 0 and 1 specifying the sparsity level in the PLS regression. Default is 1 (no sparsity).
#' @param rmax An integer specifying the maximum number of factors to consider when \code{ncomp} is not defined. Default is 5.
#' @param seed An integer seed for reproducibility. Default is 1.
#'
#' @return A matrix of generated knockoff variables of \eqn{n \times p}
#'
#' @examples
#' # Example 1: Default usage
#' library(PLSKO)
#' set.seed(123)
#' n = 100
#' p = 10
#' X <- matrix(rnorm(n*p), n, p)
#' Xk <- plsko(X)
#'
#' # Example 2: User provided neighbourhood list
#' set.seed(123)
#' n = 100
#' p = 10
#' X <- matrix(rnorm(n*p), n, p)
#' ## When the neighbourhood list of variables is provided as a symmetric adjacency matrix (although in this example all variables in X are independent)
#' p_nb <- 0.1 # probability of having a neighbour
#' nb.adj <- matrix(rbinom(p*p, 1, p_nb/2), p, p)
#' nb.adj <- nb.adj + t(nb.adj) # make it symmetric
#' nb.adj <- ifelse(nb.adj!=0, 1, 0) # make sure it is binary
#' diag(nb.adj) <- 0  # make sure diagonal is not a neighbour
#' Xk <- plsko(X, nb.list = nb.adj)
#'
#' ## When the neighbourhood list of variables is provided as a list
#' set.seed(123)
#' n = 100
#' p = 10
#' X <- matrix(rnorm(n*p), n, p)
#' ### Generate a random nb.list with 1-3 neighbours for each variable (although in this example all variables in X are independent)
#' nb.list <- list()
#' for(i in 1:p){
#'  nb.list[[i]] <- sample(1:p, sample(1:3, 1), replace = F)
#'  nb.list[[i]] <- nb.list[[i]][nb.list[[i]] != i] # exclude itself
#'  }
#'  Xk <- plsko(X, nb.list = nb.list)
#'
#' @references Yang G et al. PLSKO: a robust knockoff generator to control false discovery rate in omics variable selection. 2024:2024.08.06.606935.
#' @export

plsko <- function(X, nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = NULL, sparsity = 1, rmax = 5, seed = 1){
  set.seed(seed)

  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)

  #Input type validation
  # Step 1: Validate X (Design matrix)
  if (is.data.frame(X)) {
    X.names <- names(X)
    X <- as.matrix(X)
  } else if (is.matrix(X)) {
    X.names <- colnames(X)
  } else {
    stop("Input X must be a numeric matrix or data frame.")
  }
  if (!is.numeric(X)) stop("Input X must be a numeric matrix or data frame.")

  # Step 2: Validate nb.list (Neighborhood list or adjacency matrix)
  if (!is.null(nb.list)) {
    if (!(is.matrix(nb.list) || is.list(nb.list))) {
      stop("Input nb.list must be a list or a matrix.")
    }

    if (is.matrix(nb.list)) {
      # Ensure the matrix is symmetric, numeric, binary, and p x p
      if (!isSymmetric(nb.list)) stop("Input nb.list must be a symmetric matrix.")
      if (!is.numeric(nb.list)) stop("Input nb.list must be a numeric matrix.")
      if (any(nb.list != 0 & nb.list != 1)) stop("Input nb.list must be a binary matrix.")
      if (!all(dim(nb.list) == c(p, p))) stop("Input nb.list must be a p x p matrix.")
      if (any(diag(nb.list) != 0)) {
        warning("Input nb.list must have 0 on the diagonal. Setting diagonal to 0.")
        diag(nb.list) <- 0
      }
    }

    if (is.list(nb.list) && length(nb.list) != p) {
      stop("Input nb.list must be a list of length p.")
    }
  } else {
    # Step 3: Validate threshold.abs and threshold.q
    if (!is.null(threshold.abs)) {
      if (!is.numeric(threshold.abs) || threshold.abs < 0 || threshold.abs > 1) {
        stop("Input threshold.abs must be a numeric value between 0 and 1.")
      }
    }

    if (!is.null(threshold.q)) {
      if (!is.numeric(threshold.q) || threshold.q < 0 || threshold.q > 1) {
        stop("Input threshold.q must be a numeric value between 0 and 1.")
      }
    }
  }

# If nb.list is not provided, generate a neighborhood list based on correlations
  if(is.null(nb.list)){
    #random swap column of X
    sample.order.mat <- diag(1, p, p)[,sample(p)]
    X <- X %*% sample.order.mat

    mu <- colMeans(X)
    call <- list(thres.abs = threshold.abs,
                 thres.q = threshold.q,
                 ncomp = ncomp,
                 sparsity = sparsity)

    X <- scale(X, center = T, scale = F)

    # Calculate the correlation matrix and replace diagonal with 0
    cor_mat <- cor(X)
    diag(cor_mat) <- 0

    # Determine the correlation threshold
    if(!is.null(threshold.abs)){
      threshold <- threshold.abs
    }
    else if(!is.null(threshold.q)){
      threshold <- quantile(unlist(abs(cor_mat)), prob = threshold.q, names = F)
    }
    else {
      threshold <- quantile(unlist(abs(cor_mat)), prob = 0.80, names = F)
    }

    # Create a list of neighborhoods based on the correlation threshold
    neighborhoods <- vector(mode = "list", length = p)
    for (i in 1:p) {
      # Find the indices of the columns with correlation greater than the threshold
      cols <- which(abs(cor_mat[i,]) > threshold)

      # Add the indices to the list of neighborhoods
      if (length(cols) > 0) {
        neighborhoods[[i]] <- cols
      }
    }
  }
  else{
    # If a neighborhood list is provided by the user
    if("list" %in% class(nb.list)){
      neighborhoods <- nb.list
    }
    else if("matrix" %in% class(nb.list) & isSymmetric(nb.list)){
      neighborhoods <- vector(mode = "list", length = p)
      for (i in 1:p) {
        # Find the indices of the columns with correlation greater than the threshold
        cols <- which(nb.list[i,] == 1)
        # Add the indices to the list of neighborhoods
        if (length(cols) > 0) {
          neighborhoods[[i]] <- cols
        }
      }
    }
  }

  # Generate knockoff variables using PLS regression
  ## if neighbours more than or equal to 2, using PLS regression to estimate variable X_i's conditional distribution;
  ## if neighbours less than 2, using linear regression
  ## if no neighbours, permute X_i
  X <- as.data.frame(X)
  X_k <- matrix(NA,nrow(X),ncol(X))
  X_k <- as.data.frame(X_k)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(colnames(X),"k")

  # If ncomp is not provided, set it the minimum of p/2 and the empirical number of components
  if(is.null(ncomp)){
    r_emp <- r_criterion(X, rmax = rmax)
  }

  # Initialize the progress bar
  pb <- progress_bar$new(
    format = " Generating knockoff variables [:bar] :percent in :elapsed",
    total = p, clear = FALSE, width = 60
  )

  # Loop over each variable to estimate its conditional distribution
  for (i in 1:p){
    nb <- neighborhoods[[i]][neighborhoods[[i]]!=i] #exclude itself!
    nb <- nb[!duplicated(nb)] # remove duplicates
    k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
    X.nb <- X[,nb]
    X_k.nb <- X_k[, k.nb]

    if(length(nb) ==0){
      Y <- X[,i]
      Y.hat <- 0
    }

    else if(length(nb)== 1){
      X.run <- cbind(X.nb, X_k.nb)
      Y <- X[,i]

      Y.hat <- linear.regression.generator(Y, X.run)
    }

    else{
      X.run <- cbind(X.nb, X_k.nb)
      Y <- X[,i]
      if(is.null(ncomp)){ #if ncomp is not provided, set it the minimum of p/2 and the empirical number of components
        this.ncomp <- ceiling(min(ncol(X.nb)/2, r_emp))
      }
      else this.ncomp <- ncomp

      this.ncomp <- max(this.ncomp, 2) #minimum 2 components

      #when sparsity < 1, sparse PLS regression is used for conditional distribution with sparse*p kept on each comp
      keepX <- rep(round(sparsity*ncol(X.run)), this.ncomp)
      Y.hat <- pls.recovery.generator(Y, X.run, ncomp = this.ncomp, keepX = keepX)
    }

    # calculate the residuals and permute
    Y.res <- Y - Y.hat
    res <- sample(Y.res)

    X_k[,i] <- Y.hat + res

    # Update the progress bar
    pb$tick()

  }

  # Add the mean back to knockoff variables generated from the centered data
  X_k <- apply(X_k, 1, function(x){x+mu})
  X_k <- t(X_k)

  # Swap columns back to original order if necessary
  if(is.null(nb.list)){
    #swap back
    X_k <- X_k %*% t(sample.order.mat)
  }

  #obj <- list(X_k = X_k,
  #            call = call)
  return(X_k)
}


#' Calculate \eqn{\hat{X}} by fitting PLS regression on its neighbours
#'
#' @rdname pls.recovery.generator
#' @keywords internal
#'
pls.recovery.generator <- function(Y, X, ncomp, keepX = rep(ncol(X), ncomp)){
  X <- as.data.frame(X)
  #X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- mixOmics::spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")
  Y.hat <- predict(pls, X)$predict[1:n, 1, ncomp]
  return(Y.hat)
}

#' Calculate \eqn{\hat{X}} by fitting OLS regression on its neighbours
#'
#' @rdname linear.regression.generator
#' @keywords internal
#'
linear.regression.generator <- function(Y, X){
  X <- as.data.frame(X)
  X <- X[!duplicated(as.list(X))]
  X <- as.matrix(X)
  betas <- solve(crossprod(X))%*% t(X) %*% Y
  Y.hat <- X %*% betas
  return(Y.hat)
}


