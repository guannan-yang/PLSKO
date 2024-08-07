library(mixOmics)
plsko <- function(X, nb.list = NULL, threshold.abs = NULL, threshold.q = NULL, ncomp = 2, sparsity = 1){
  n <- nrow(X)
  p <- ncol(X)
 
  if(is.null(nb.list)){
    
    #random swap column of X
    sample.order.mat <- diag(1, p, p)[,sample(p)]
    X <- X %*% sample.order.mat
    
    mu <- colMeans(X)  
    call <- list(thres.abs = threshold.abs,
                 thres.q = threshold.q,
                 ncomp = ncomp,
                 sparsity = sparsity,
                 gauss.err = gauss.err)
    
    X <- scale(X, center = T, scale = F)
    
    # Calculate the correlation matrix
    cor_mat <- cor(X)
    # Replace the diagonal values with 0, since they are always 1
    diag(cor_mat) <- 0
    
    if(!is.null(threshold.abs)){
      threshold <- threshold.abs
    }  
    else if(!is.null(threshold.q)){
      threshold <- quantile(unlist(abs(cor_mat)), prob = threshold.q, names = F)
    } 
    else {
      threshold <- quantile(unlist(abs(cor_mat)), prob = 0.90, names = F)
    }
    
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
  else{ #user provided neighbour list
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
        this.ncomp <- ceiling(min(n/2, ncol(X.nb)/2))
      }
      else this.ncomp <- ncomp
      
      this.ncomp <- min(this.ncomp, ncol(X.run)-1)
      
      #when sparsity < 1, sparse PLS regression is used for conditional distribution with sparse*p kept on each comp
      keepX <- rep(round(sparsity*ncol(X.run)), this.ncomp) 
      Y.hat <- pls.recovery.generator(Y, X.run, ncomp = this.ncomp, keepX = keepX)
    }
    
    # calculate the residuals
    Y.res <- Y - Y.hat
    
    # permuting the residuals
    res <- sample(Y.res) 
      
    X_k[,i] <- Y.hat + res
  }
  
  X_k <- apply(X_k, 1, function(x){x+mu})
  X_k <- t(X_k)
  
  if(is.null(nb.list)){
    #swap back
    X_k <- X_k %*% t(sample.order.mat)
  }

  obj <- list(X_k = X_k, 
              call = call)
  return(obj)  
}


pls.recovery.generator <- function(Y, X, ncomp, keepX = rep(ncol(X))){
  X <- as.data.frame(X)
  #X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- spls(X, Y, ncomp = ncomp, scale = F, keepX = keepX, mode = "regression")
  Y.hat <- predict(pls, X)$predict[1:n, 1, ncomp]
  return(Y.hat)
}

linear.regression.generator <- function(Y, X){
  X <- as.data.frame(X)
  X <- X[!duplicated(as.list(X))]
  X <- as.matrix(X)
  betas <- solve(crossprod(X))%*% t(X) %*% Y
  Y.hat <- X %*% betas
  return(Y.hat)
}

# calculate the empirical number of component used for PLS regression by the PC_p1 criterion from Bai and Ng (2002).
r_criterion <- function(X, rmax = 20){
  # rmax = 10 # maximum number of factors
  n = nrow(X)
  p = ncol(X)
  
  PC = matrix(NA,rmax+1,1)
  SX = scale(X)
  mXX = SX %*%t(SX)
  for(k in rmax:0){
    if (k == 0){
      PC[k+1,] = sum(SX^2/(n*p))
    } else 
      eig <- RSpectra::eigs_sym(mXX, k, which = "LM", sigma = NULL, lower = TRUE)
    meigvec <- eig$vectors 
    mF <- sqrt(n) * meigvec       # estimated factors
    Lam <- (t(mF) %*% SX)/n
    if(k==rmax){
      sigma2 = sum((SX-mF %*% Lam)^2)/(n*p)
    }
    PC[k+1,] = sum((SX-mF %*% Lam)^2)/(n*p) + k*sigma2*((n+p)/(n*p))*log((n*p)/(n+p))
  }
  r = which(PC == min(PC))[1]-1
}
