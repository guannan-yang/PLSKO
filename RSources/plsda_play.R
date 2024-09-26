# Playing around 
set.seed(7749)
plsda.for.fun <- function(X, threshold = NULL){
  n <- nrow(X)
  p <- ncol(X)
  mu <- colMeans(X)  
  
  X <- scale(X, center = T, scale = F)

  # Calculate the correlation matrix
  cor_mat <- cor(X)
  # Replace the diagonal values with 0, since they are always 1
  diag(cor_mat) <- 0
  
  if(is.null(threshold)){
    threshold <- quantile(unlist(abs(cor_mat)), prob = 0.95, names = F)
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

  # if neighbours more than or equal to n, using plsda to estimate variable i's distribution;
  # if neighbours less than n, using linear regression
  # if no neighbours, treat as random gaussian then sampling
  X <- as.data.frame(X)
  X_k <- matrix(NA,nrow(X),ncol(X))
  X_k <- as.data.frame(X_k)
  rownames(X_k) <- rownames(X)
  colnames(X_k) <- paste0(colnames(X),"k")
 for (i in 1:p){
    tryCatch({  
    print(i)
      nb <- neighborhoods[[i]]
      k.nb <- neighborhoods[[i]][neighborhoods[[i]] < i]
      X.nb <- X[,nb]
      X_k.nb <- X_k[, k.nb]
      
      if(length(nb) ==0){
        print("0")
        Y <- X[,1]
        Y.hat <- 0
      }
      
      else if(length(nb)< n/2){
        print("linear")
        X.run <- cbind(X.nb, X_k.nb)
        Y <- X[,i]

        Y.hat <- linear.regression.generator(Y, X.run)
      }
      
      else{
        if(length(k.nb) < n/2){
          print("pls")
          X.run <- cbind(X.nb, X_k.nb)
          Y <- X[,i]
          ncomp <- round(n/2)-1 # avoiding overfitting
          Y.hat <- plsda.recovery.generator(Y, X.run, ncomp = 2)
        }
        else{
          print("block.pls")
          X.run <- list(X.nb = X.nb, X_k.nb = X_k.nb)
          Y <- cbind(X[,i],X[,i])
          Y <- as.matrix(Y)
          rownames(Y) <- rownames(X)
          
          #block.pls require Y to be a matrix. won't affect results
          ncomp <- round(min(length(k.nb), n)/2)
          
          #design matrix
          Y.weight <- c(length(nb)/(length(nb)+length(k.nb)), length(k.nb)/(length(nb)+length(k.nb)))
          des.mat <- matrix(c(0,1,Y.weight[1],1,0,Y.weight[2],Y.weight[1],Y.weight[2],0),3,3)
          colnames(des.mat) <- c("X.nb","X_k.nb","Y")
          rownames(des.mat) <- c("X.nb","X_k.nb", "Y")
          
          Y.hat <- block.pls.recovery.generator(Y, X.run, ncomp = 2, design = des.mat)
          Y <- X[,i] # recovery to vector
        }
        
      }
      
      Y.res <- Y - Y.hat
      # resSD <- sqrt(sum(Y.res^2))
      # res <- rnorm(n, sd = resSD)
      # X_ki <- Y.hat + resSD
      # permuation
      res.perm <- sample(Y.res)
      X_k[,i] <- Y.hat + res.perm
    }
    , error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
   }
  
  X_k <- apply(X_k, 1, function(x){x+mu})
  X_k <- t(X_k)
  
  return(X_k)  
}

block.pls.recovery.generator <- function(Y, X, ncomp, design){
  n = nrow(X[[1]])
  pls <- block.pls(X, Y, ncomp = ncomp, scale = F, mode = "regression")
  Y.hat <- predict(pls, X)[["WeightedPredict"]][1:n, 1, ncomp]

  return(Y.hat)
}

plsda.recovery.generator <- function(Y, X, ncomp){
  X <- as.data.frame(X)
  X <- X[!duplicated(as.list(X))]
  n <- nrow(X)
  p <- ncol(X)
  pls <- pls(X, Y, ncomp = ncomp, scale = F, mode = "regression")
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


# ko.iter$PLSKO$default <- list()
# set.seed(7749)
# for(i in 1:10){
#   ko.iter$PLSKO$default[[i]] <- cov.ko(dt, method1 = "PLSKO", method2 = "default")
# }