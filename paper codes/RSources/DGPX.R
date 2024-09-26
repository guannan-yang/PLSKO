#source("/data/gpfs/projects/punim0613/Guan/Sim2/RSource/data_sim1.R")
#setwd("~/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/PhD/vitual_desk/Sim2")
source("./RSource/data_sim1.R")
DGPX.gauss <- function(n, p.in.cluster, p.indep = 0, cluster = 1, 
                       X.method = "AR1", rho = 0.5, r.in.cluster = 3, load = c(1,2), noise = 0,
                       seed = Sys.time()){
  
  set.seed(seed)
  p.total = p.in.cluster * cluster + p.indep
  r = r.in.cluster
  
  X.full <- NULL
  cov.full <- matrix(0, p.total, p.total)
  diag(cov.full) <- 1
  index.in.this.cluster <- 1:p.in.cluster
  
  if(length(X.method)==1) X.method <- rep(X.method, cluster)
  
  for(i.cluster in 1:cluster){
    method.i <- X.method[i.cluster]
    
    if(method.i == "fac1"){
      Fa <- matrix(rnorm(n*r),n,r)
      Lam <- matrix(rnorm(r*p.in.cluster),r,p.in.cluster)
      X <- Fa %*% Lam
      cov <- t(Lam) %*% Lam
    }
    
    else if(method.i == "fac2"){
      Fa <- matrix(rnorm(n*r),n,r)
      Lam.pool <- c(runif(r*p.in.cluster, min = load[1], max = load[2]), runif(r*p.in.cluster, min = -load[2], max = -load[1]))
      Lam <- matrix(sample(Lam.pool, r*p.in.cluster, replace = F), r, p.in.cluster)
      X <- Fa %*% Lam
      cov <- t(Lam) %*% Lam
    }
    else if(method.i == "AR1"){
      obj <- norm.sim(n, p.in.cluster, rho = rho)
      X <- obj$X
      cov <- obj$cov
    }
    
    else if(method.i == "equi"){
      obj <- equi.sim(n, p.in.cluster, rho = rho)
      X <- obj$X
      cov <- obj$cov
    }
    
    # else if(method.i == "sparse"){
    #   obj <- sparse.gauss.sim
    # }
    
    else if(method.i == "quad"){
      p.half <- round(p.in.cluster/2)
      X_1 <- MASS::mvrnorm(n, mu = rep(0, p.half), Sigma = diag(1, nrow = p.half, ncol = p.half))
      X_2 <- X_1^2
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    else if(method.i == "fac2-quad"){
      p.half <- round(p.in.cluster/2)
      Fa <- matrix(rnorm(n*r),n,r)
      Lam.pool <- c(runif(r*p.half, min = load[1], max = load[2]), runif(r*p.half, min = -load[2], max = -load[1]))
      Lam <- matrix(sample(Lam.pool, r*p.half, replace = F), r, p.half)
      X_1 <- Fa %*% Lam
      X_2 <- X_1^2
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    else if(method.i == "equi-quad"){
      p.half <- round(p.in.cluster/2)
      obj <- equi.sim(n, p.half, rho = rho)
      X_1 <- obj$X
      X_2 <- X_1^2
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    else if(method.i == "fac2-sin"){
      p.half <- round(p.in.cluster/2)
      Fa <- matrix(rnorm(n*r),n,r)
      Lam.pool <- c(runif(r*p.half, min = load[1], max = load[2]), runif(r*p.half, min = -load[2], max = -load[1]))
      Lam <- matrix(sample(Lam.pool, r*p.half, replace = F), r, p.half)
      X_1 <- Fa %*% Lam
      X_2 <- sin(2*pi/max(X_1)*X_1)
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    else if(method.i == "fac2-inter"){
      p.half <- round(p.in.cluster/2)
      Fa <- matrix(rnorm(n*r),n,r)
      Lam.pool <- c(runif(r*p.half, min = load[1], max = load[2]), runif(r*p.half, min = -load[2], max = -load[1]))
      Lam <- matrix(sample(Lam.pool, r*p.half, replace = F), r, p.half)
      X_1 <- Fa %*% Lam
      # Generate the interaction variables
      X_2 <- matrix(0, nrow = n, ncol = p.half)
      for (i in 1:(p.half - 1)) {
        X_2[, i] <- X_1[, i] * X_1[, i + 1]
      }
      # Generate the last interaction variable as the product of the last and first variable in the independent half
      X_2[, p.half] <- X_1[, p.half] * X_1[, 1]
      
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    else if(method.i == "equi-inter"){
      p.half <- round(p.in.cluster/2)
      obj <- equi.sim(n, p.half, rho = rho)
      X_1 <- obj$X
      # Generate the interaction variables
      X_2 <- matrix(0, nrow = n, ncol = p.half)
      for (i in 1:(p.half - 1)) {
        X_2[, i] <- X_1[, i] * X_1[, i + 1]
      }
      # Generate the last interaction variable as the product of the last and first variable in the independent half
      X_2[, p.half] <- X_1[, p.half] * X_1[, 1]
      
      X <- cbind(X_1, X_2)
      cov <- NA
    }
    
    X.full <- cbind(X.full, X)
    cov.full[index.in.this.cluster, index.in.this.cluster] <- cov
    
    index.in.this.cluster <- index.in.this.cluster + p.in.cluster
  }

  if(p.indep !=0 ){
  X.indep <- MASS::mvrnorm(n, rep(0, p.indep), Sigma = diag(nrow = p.indep))
  X.full <- cbind(X.full, X.indep)
  cov.full[index.in.this.cluster + p.indep, index.in.this.cluster+p.indep] <- diag(x = 1, nrow = p.indep)
  }

  X.full <- X.full + matrix(rnorm(n*p.total, 0, sd = noise), nrow = n, ncol = p.total)
  cov.full <- cov.full + diag(x = noise^2, nrow = p.total, ncol = p.total)
  
  return(list(X = X.full,
              cov = cov.full))
}

  
equi.sim <- function(n, p, rho){
  Sigma <- matrix(rho,nrow = p, ncol = p)
  diag(Sigma) <- 1
  
  mu <- rep(0, p)
  X <- MASS::mvrnorm(n, mu=mu, Sigma=Sigma)
  u <- pnorm(X)
  
  X.norm <- qnorm(u)
  return(list(X = X.norm,
              cov = Sigma))  
  
}
  