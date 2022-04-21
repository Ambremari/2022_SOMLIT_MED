##REGRESSION POLYNOMIALE LOCALE

reg_pol <- function(Xech, yech, h, p, k=1e-2){
  Xech <- Xech*k
  h <- h*k
  n <- length(Xech)
  ti <- min(Xech)
  tf <- max(Xech)
  X <- seq(ti, tf, k)
  m <- length(X)
  Mp <- NULL
  for(i in 1:m){ 
    Xx <- NULL
    for(o in 0:p){
      xx <- matrix((X[i]-Xech)^o, ncol=1)
      Xx <- cbind(Xx, xx)
    }
    Wx <- poids_x(X[i], Xech, h)
    Wx <- 1/h*diag(Wx, nrow=n)
    A <- t(Xx)%*%Wx%*%Xx
    beta_hat <- solve(A) %*% t(Xx) %*% Wx %*% yech
    e1 <- matrix(c(1, rep(0,p)), ncol=1)
    Sx <- t(e1)%*%solve(A) %*% t(Xx) %*% Wx 
    beta0_hat <- beta_hat[1,1]
    Mp <- c(Mp, beta0_hat)
  }
  return(list(X=X*k^(-1), Mp=Mp))
}
