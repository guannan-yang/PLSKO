DGPX.fac <- function(n, p.in.cluster, p.indep = 0, cluster = 1, 
                     X.method = "fac2", rho = 0.5, r.in.cluster = 3, load = c(1,2), noise = 0,
                     seed = Sys.time()){
  
  set.seed(seed)
  p.total = p.in.cluster * cluster + p.indep
  r = r.in.cluster
  
  X.full <- NULL
  cov.full <- matrix(0, p.total, p.total)
  diag(cov.full) <- 1
  
  Lam.full <- matrix(0, r*cluster, p.total)
  Fa.full <- matrix(0, n, r*cluster)
  
  index.in.this.cluster <- 1:p.in.cluster
  rank.in.this.cluster <- 1:r
  
  if(length(X.method)==1) X.method <- rep(X.method, cluster)
  
  for(i.cluster in 1:cluster){
    method.i <- X.method[i.cluster]
    
    if(method.i == "fac1"){
      Fa <- matrix(rnorm(n*r),n,r)
      Fa <- Fa/sqrt(n-1)
      Lam <- matrix(rnorm(r*p.in.cluster),r,p.in.cluster)
      X <- Fa %*% Lam
      cov <- t(Lam) %*% Lam
    }
    
    else if(method.i == "fac2"){
      Fa <- matrix(rnorm(n*r),n,r)
      Fa <- Fa/sqrt(n-1)
      Lam.pool <- c(runif(r*p.in.cluster, min = load[1], max = load[2]), runif(r*p.in.cluster, min = -load[2], max = -load[1]))
      Lam <- matrix(sample(Lam.pool, r*p.in.cluster, replace = F), r, p.in.cluster)
      X <- Fa %*% Lam
      cov <- t(Lam) %*% Lam
    }
   
    X.full <- cbind(X.full, X)
    cov.full[index.in.this.cluster, index.in.this.cluster] <- cov

    Lam.full[rank.in.this.cluster, index.in.this.cluster] <- Lam
    Fa.full[, rank.in.this.cluster] <- Fa
    
    index.in.this.cluster <- index.in.this.cluster + p.in.cluster
    rank.in.this.cluster <- rank.in.this.cluster + r
  }
  
  if(p.indep !=0 ){
    X.indep <- MASS::mvrnorm(n, rep(0, p.indep), Sigma = diag(nrow = p.indep))
    X.full <- cbind(X.full, X.indep)
    cov.full[index.in.this.cluster + p.indep, index.in.this.cluster+p.indep] <- diag(x = 1, nrow = p.indep)
  }
  
  X.full <- X.full + matrix(rnorm(n*p.total, 0, sd = noise), nrow = n, ncol = p.total)
  cov.full <- cov.full + diag(x = noise^2, nrow = p.total, ncol = p.total)
  
  return(list(X = X.full,
              cov = cov.full,
              Lam = Lam.full,
              Fa = Fa.full))
}

DGPy.fac <- function(X.list, s, A, c = 0, y.dis = "Normal"){
  n = nrow(X)
  p = ncol(X)
  B.true <- rep(0,p)
  S.true <- sample(1:p, size = s) 
  S.true <- sort(S.true) # True Support
  C.unif <- c(runif(p, A[1], A[2]), runif(p, -A[2], -A[1]))
  C <- sample(C.unif, size = s , replace = TRUE) # coefficients for true signals
  B.true [S.true] <- C # A true vector of coefficients
  epsilon <- matrix(rnorm(n),n,1)  # Error Term
  
  linear.fit <- X %*% B.true # Response variable
  
  if(y.dis == "Normal"){
    y <- linear.fit + sqrt(c*s)*epsilon
    # noise.ratio = var(sqrt(c*s)*epsilon)/var(linear.fit)
    # y = list(y, noise.ratio)
  }
  else if(y.dis == "NB"){
    mu = 2^linear.fit*400
    prob = 0.97
    size = prob/(1-prob)*mu
    y <- rnbinom(n, size = size, mu = mu)
  }
  
  else if(y.dis == "Binary"){    
    p.logit <- 1/(1+exp(-(linear.fit-mean(linear.fit))))# logistics
    y <- rbinom(n, size = 1, prob = p.logit)
  }
  
  obj = list(Beta = B.true, y = y)
  # }
  return(obj)
}