#' @title Generate multiple PLSKO knockoff variables
#' @description This function generates multiple knockoff variables for a given design matrix using the PLS regression method with mutliple times residuals only. Please note this function is fast but does not generate rigidly valid knockoff variables under the definition.
#' @param X A numeric matrix or data frame where rows represent observations and columns represent variables.
#' @param n_ko An integer specifying the number of knockoff variables to generate. Default is 10.
#'
#' @return A list of length n_ko, each element is a matrix of n*p of knockoff variables.
#'
#' @keywords internal

plsko_batch <- function(X, n_ko = 10, nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = NULL, sparsity = 1, rmax = 5, seed = Sys.time(), simpls = F){
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

  # genearte X_k.list of length n_ko, each element is a matrix of n*p of NA
  X_k.list <- vector(mode = "list", length = n_ko)
  for(j in 1:n_ko){
    X_k.list[[j]] <- matrix(NA,nrow(X),ncol(X))
    X_k.list[[j]] <- as.data.frame(X_k.list[[j]])
    rownames(X_k.list[[j]]) <- rownames(X)
    colnames(X_k.list[[j]]) <- paste0(colnames(X),"k")
  }

  # Initialize the progress bar
  pb <- progress::progress_bar$new(
    format = " Generating knockoff variables [:bar] :percent in :elapsed",
    total = p, clear = FALSE, width = 60
  )

  # Loop over each variable to estimate its conditional distribution
  for (i in 1:p){
    nb <- neighborhoods[[i]][neighborhoods[[i]]!=i] #exclude itself!
    nb <- nb[!duplicated(nb)] # remove duplicates
    k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
    X.nb <- X[,nb]

    if(length(nb) ==0){
      Y <- X[,i]
      Y.hat <- 0
    }

    else if(length(nb)== 1){
      X.run <- X.nb
      Y <- X[,i]

      Y.hat <- linear.regression.generator(Y, X.run)
    }

    else{
      X.run <- X.nb
      Y <- X[,i]
      if(is.null(ncomp)){ #if ncomp is not provided, set it the minimum of p/2 and the empirical number of components
        r_emp <- r_criterion(X, rmax = rmax)
        this.ncomp <- ceiling(min(ncol(X.nb)/2, r_emp))
      }
      else{
        this.ncomp <- ncomp
      }

      this.ncomp <- min(this.ncomp, ncol(X.run)) #maximum ncomp is the number of variables minus one in the regression
      this.ncomp <- max(this.ncomp, 2) #minimum 2 components

      #when sparsity < 1, sparse PLS regression is used for conditional distribution with sparse*p kept on each comp
      keepX <- rep(round(sparsity*ncol(X.run)), this.ncomp)
      Y.hat <- pls.recovery.generator(Y, X.run, ncomp = this.ncomp, keepX = keepX)
    }

    # calculate the residuals and permute
    Y.res <- Y - Y.hat

    for(j in 1:n_ko){
    res <- sample(Y.res)

    X_k.list[[j]][,i] <- Y.hat + res
    }
    # Update the progress bar
    pb$tick()

  }

  # Add the mean back to knockoff variables generated from the centered data
  X_k <- apply(X_k, 1, function(x){x+mu})
  X_k <- t(X_k)
  X_k.list <- lapply(X_k.list, function(x) {a <- apply(x, 1, function(y){y+mu}); t(a)})

  # Swap columns back to original order if necessary
  if(is.null(nb.list)){
    #swap back
    X_k <- X_k %*% t(sample.order.mat)
    X_k.list <- lapply(X_k.list, function(x) {a <- x %*% t(sample.order.mat); a})
  }


  #obj <- list(X_k = X_k,
  #            call = call)
  return(X_k.list)
}
