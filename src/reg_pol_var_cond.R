##REGRESSION POLYNOMIALE LOCALE
#Variance conditionnelle

var_cond <- function(Xech, yech, h, p, a, sigma, k=1e-2){
  Xech <- Xech*k
  h <- h*k
  n <- length(Xech)
  ti <- min(Xech)
  tf <- max(Xech)
  X <- seq(ti, tf, k)
  m <- length(X)
  #var sorties
  mat_var <- NULL
  for(i in 1:m){ 
    Xx <- NULL
    for(o in 0:p){
      xx <- matrix((X[i]-Xech)^o, ncol=1)
      Xx <- cbind(Xx, xx)
    }
    Wx <- poids_x(X[i], Xech, h)
    Wx <- 1/h*diag(Wx, nrow=n)
    A <- t(Xx)%*%Wx%*%Xx
    mv <- solve(A)%*%(t(Xx)%*%Wx^2%*%Xx)%*%solve(A)
    mv <- mv*sigma[i]
    mat_var <- cbind(mat_var, mv)
  }
  mat_var <- array(mat_var, dim=c(p+1, p+1, m)) #Matrices de variance conditionelle p+1 x p+1 pour m x
  return(Mat_var=mat_var)
}
