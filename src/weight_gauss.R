#POIDS NOYAU GAUSSIEN

poids_x <- function(x, xech, h){
  u <- 1/h*(x-xech)
  Ku <- 1/sqrt(2*pi)*exp(-1/2*u^2) 
  return(Ku)
}