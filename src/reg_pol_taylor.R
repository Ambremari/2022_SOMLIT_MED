#Régression polynomiale locale
#Approximation de Taylor

reg_pol_taylor <- function(Xech, yech, h, p, a, k=1e-2){
  h <- h*k
  n <- length(Xech)
  ti <- min(Xech)
  tf <- max(Xech)
  X <- seq(ti, tf, 1)
  Xech <- Xech*k
  X <- X*k
  m <- length(X)
  #var sortie
  Bias <- NULL
  Mp <- NULL
  V <- NULL
  all_denom <- NULL
  all_Kx <- NULL
  all_dn <- NULL
  ord <- p+a
  for(i in 1:m){ 
    Xx <- NULL
    for(o in 0:ord){
      xx <- matrix((X[i]-Xech)^o, ncol=1)
      Xx <- cbind(Xx, xx)
    }
    #regression polynomiale
    Kx <- poids_x(X[i], Xech, h)
    Wx <- 1/h*diag(Kx, nrow=n)
    A <- t(Xx)%*%Wx%*%Xx
    mv <- solve(A)%*%(t(Xx)%*%Wx^2%*%Xx)%*%solve(A) 
    beta_hat <- solve(A) %*% t(Xx) %*% Wx %*% yech
    beta0_hat <- beta_hat[1,1]
    Mp <- c(Mp, beta0_hat) #estimation
    #V
    mv <- solve(A)%*%(t(Xx)%*%Wx^2%*%Xx)%*%solve(A)
    V <- cbind(V, mv)
    #degrés de liberté
    dn <- sum(Kx)^2/sum(Kx^2)
    all_dn <- c(all_dn, dn)
    #biais Sn-1
    tau <- seq(p+1, 2*p+a, 1)
    Snj <- NULL
    for(j in tau){
      J <- matrix((X[i]-Xech)^j, ncol=1)
      SN <- Kx %*% J
      Sn <- ifelse(j<=ord, sum(SN), 0)
      Snj <- cbind(Snj, Sn)
    }
    betaa <- matrix(beta_hat[(p+2):(ord+1)], ncol=1)
    mat_Sn <- NULL
    for(l in 1:(1+p)){
    row_l <- Snj[l:(l+a-1)]
    mat_Sn <- rbind(mat_Sn, row_l)
    }
    bias <- mat_Sn%*%betaa
    Bias <- cbind(Bias, bias) #biais
    #somme des carrées des résidus pondérés normalisée
    denom <- Wx - Wx%*%Xx%*%solve(A)%*%t(Xx)%*%Wx
    denom <- sum(diag(denom))
    all_denom <- c(all_denom, denom)
    all_Kx <- cbind(all_Kx, Kx)
  }
  V <- array(V, dim=c(ord+1, ord+1, m))
  yind <- which(X %in% Xech)
  yestim <- Mp[yind]
  res <- (yech-yestim)^2
  Sigma <- NULL
  for(i in 1:m){
    nume <- res %*% all_Kx[,i]
    nume <- sum(nume)
    sigma <- nume/all_denom[i]
    Sigma <- c(Sigma, sigma)
  }
  return(list(X=X*k^(-1), Mp=Mp, biais=Bias, sigma=Sigma, V=V, dn=all_dn))
}
