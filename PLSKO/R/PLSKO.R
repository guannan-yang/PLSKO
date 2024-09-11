#' Generate a knockoff variable set with PLSKO using PLS regression
#'
#' @importFrom mixOmics spls
#' @param X A numeric matrix or data frame. The original design data matrix with \eqn{n} observations as rows and \eqn{p} variables as columns.
#' @param nb.list Optional. A list of length \eqn{p} or adjacency matrix of \eqn{p \times p} that defines the neighbourship of variables. A list of length \eqn{p} should include the neighbours' index of each variable from \eqn{X_1} to \eqn{X_p} in order; The \eqn{i^{th}} element in the list includes the indices of the neighbour variables of \eqn{X_i}, or \code{NULL} when no neighbours. A adjacency matrix should be symmetric with only binary element and  where \eqn{M_{ij} = 1} when \eqn{X_i} and \eqn{X_j} is defined as neighbours; otherwise \eqn{M_{ij} = 0} when not neighbour or on diagonal (i.e. \eqn{i = j}). If not provided or NULL, the neighborhoods are determined based on correlations.
#' @param threshold.abs Optional. A value between \eqn{0} and \eqn{1}. A numeric value specifying an absolute correlation threshold to define neighborhoods. If not provided, the function uses \code{threshold.q}.
#' @param threshold.q Optional. A numeric value between 0 and 1 indicating the quantile of the correlation values to use as a threshold.
#' @param ncomp Optional. An integer specifying the number of components to use in the PLS regression. Default is 2.
#' @param sparsity Optional. A numeric value between 0 and 1 specifying the sparsity level in the PLS regression. Default is 1 (no sparsity).
#'
#' @return A matrix of generated knockoff variables of \eqn{n \times p}
#' @export
#'
#' @examples
#'
plsko <- function(X, nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = NULL, sparsity = 1){
  n <- nrow(X)
  p <- ncol(X)

  #Input type validation
  if(is.data.frame(X)){
    X.name = names(X)
    X = as.matrix(X)
  }else if (is.matrix(X)) {
    X.names = colnames(X)
  }else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')

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
      threshold <- quantile(unlist(abs(cor_mat)), prob = 0.90, names = F)
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
  # if neighbours more than or equal to 2, using PLS regression to estimate variable X_i's conditional distribution;
  # if neighbours less than 2, using linear regression
  # if no neighbours, permute X_i
  X <- as.data.frame(X)
  X_k <- matrix(NA,nrow(X),ncol(X))
  X_k <- as.data.frame(X_k)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(colnames(X),"k")

  # If ncomp is not provided, set it the minimum of p/2 and the empirical number of components
  if(is.null(ncomp)){
    r_emp <- r_criterion(X)
  }

  # Loop over each variable to estimate its conditional distribution
  for (i in 1:p){
    nb <- neighborhoods[[i]]
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
      if(is.null(ncomp)){
        ncomp <- ceiling(min(ncol(X.nb)/2, r_emp))
      }
      else this.ncomp <- ncomp

      this.ncomp <- min(this.ncomp, ncol(X.run)-1)

      #when sparsity < 1, sparse PLS regression is used for conditional distribution with sparse*p kept on each comp
      keepX <- rep(round(sparsity*ncol(X.run)), this.ncomp)
      Y.hat <- pls.recovery.generator(Y, X.run, ncomp = this.ncomp, keepX = keepX)
    }

    # calculate the residuals and permute
    Y.res <- Y - Y.hat
    res <- sample(Y.res)

    X_k[,i] <- Y.hat + res
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
pls.recovery.generator <- function(Y, X, ncomp, keepX = rep(ncol(X))){
  X <- as.data.frame(X)
  #X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")
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


