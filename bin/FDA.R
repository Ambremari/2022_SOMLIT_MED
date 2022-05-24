#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#FDA
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lfe)
library(missMDA)
library(FactoMineR)
library(ggforce)

###source code###
source("src/my_palette.R")

###load data###
data_ab <- read.csv('results/PC_AB.csv')
weight_ab <- read.csv("results/PC_AB_WEIGHTS.csv")
data_nutri <- read.csv('results/PC_NUTRI.csv')
weight_nutri <- read.csv("results/PC_NUTRI_WEIGHTS.csv")
data_diff <- read.csv('results/PC_DIFF.csv')
weight_diff <- read.csv("results/PC_DIFF_WEIGHTS.csv")

###fda plots function###
my_fda <- function(data, weights, N, variable, nutri=FALSE){
  data <- data[,-1]
  data <- data[, 1:(N+1)]
  X <- data[,-1]
  ind.w <- order(weights$VAR)
  weights <- weights$WEIGHT[ind.w]
  weights <- weights[1:N]
  ##NA treatment
  imputed <- imputePCA(X, ncp=2)
  input <- imputed$completeObs
  #fda
  data <- data.frame("SITE"=data[,1], input)
  X <- input
  X <- t(t(X)*sqrt(weights))
  n <- nrow(X)
  p <- ncol(X)
  group <- data[,1]
  group <- as.factor(group)
  my_levels <- levels(group)
  counts <- as.vector(table(group))
  prop <- counts/n
  ng <- length(prop)
  group_means <- data %>% group_by(SITE) %>% summarise_all(mean)
  group_means <- as.matrix(group_means[,-1])
  gp_centered <- demeanlist(X, list(group))
  gp_centered <- as.matrix(gp_centered)
  f1 <- sqrt(diag(var(gp_centered))) #var-cov
  B.scaling <- diag(1/f1, ncol=p) #W^-(1/2)
  fac <- 1/(n-ng)
  B <- sqrt(fac) * gp_centered %*% B.scaling #B
  B.s <- svd(B)
  tol <- 1.0e-4
  rank <- sum(B.s$d > tol)
  scaling <- B.scaling %*% B.s$v[, 1L:rank] %*% diag(1/B.s$d[1L:rank], ncol=rank)
  xbar <- colSums(prop %*% group_means)
  fac <- 1/(ng-1)
  A <- sqrt((n*prop)*fac) * scale(group_means, center=xbar, scale=FALSE) %*% scaling
  A.s <- svd(A)
  rank <- sum(A.s$d>tol * A.s$d[1L])
  V <- A.s$v[,1L:rank]/sqrt(weights)
  coord.var <- scaling %*% V
  ind <- scale(X, center=xbar, scale=FALSE) %*% coord.var 
  coord.ind <- ind[, 1:rank]
  dimnames(coord.var) <- list(colnames(X), paste("LD", 1L:rank, sep=''))
  colnames(coord.ind) <- paste("LD", 1L:rank, sep='')
  contrib.var <- t(t(coord.var^2)/A.s$d[1:rank]^2)*weights
  eig <- A.s$d[1:rank]
  
  p_axe1 <- round(prop.table(eig^2)[1]*100, 2)
  p_axe2 <- round(prop.table(eig^2)[2]*100, 2)
  my_df <- data.frame("SITE"=data[,1], coord.ind)
  Xacp1 <- my_df %>% filter(SITE=="Banyuls")
  Xacp2 <- my_df %>% filter(SITE=="Marseille")
  Xacp3 <- my_df %>% filter(SITE=="Villefranche")
  g1 <- apply(Xacp1[,2:3], 2, mean)
  g2 <- apply(Xacp2[,2:3], 2, mean)
  g3 <- apply(Xacp3[,2:3], 2, mean)
  grav_center <- rbind(g1, g2, g3)
  grav_center <- data.frame("SITE"=c('Banyuls', 'Marseille', 'Villefranche'),
                            grav_center)
  ban <- my_df %>% filter(SITE=='Banyuls')
  mar <- my_df %>% filter(SITE=='Marseille')
  vil <- my_df %>% filter(SITE=='Villefranche')
  acp_plot <- my_df %>% ggplot() + geom_point(aes(LD1, LD2, col=factor(SITE))) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    theme_light() +
    geom_point(data=grav_center, aes(LD1, LD2, col=factor(SITE)), shape=3) +
    geom_text(data=grav_center, aes(LD1, LD2, col=factor(SITE), label=SITE),
              fontface='bold', size=6, vjust=1, hjust=1) +
    geom_segment(data=ban, aes(g1[1], g1[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    geom_segment(data=mar, aes(g2[1], g2[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    geom_segment(data=vil, aes(g3[1], g3[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    scale_color_manual(values=my_palette[2:4], name='Site') +
    xlab(paste('Axe 1 ', p_axe1, '%')) + ylab(paste('Axe 2 ', p_axe2, '%')) + ggtitle(variable)+
    theme(legend.text = element_text(size=16),
          legend.position='none',
          legend.title = element_text(size=18),
          axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          title=element_text(size=18))
  coord <- data.frame('VAR'=row.names(coord.var), coord.var)
  coord <- coord %>% separate('VAR', c('PC', 'GROUPE'), '_')
  leg_name <- 'Groupe'
  if(nutri==TRUE){
   leg_name <- 'Nutriment'
  }
  var_plot <- coord %>% ggplot() + 
    geom_segment(aes(0, 0, xend=LD1, yend=LD2)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), inherit.aes = FALSE, 
                size=.2, linetype='dashed', alpha=.6) +
    geom_label(aes(LD1, LD2, label=PC, fill=GROUPE),
               size=4.5, fontface='bold', alpha=.8, col='white') +
    theme_light() + xlab(paste('Axe 1 ', p_axe1, '%')) + ylab(paste('Axe 2 ', p_axe2, '%')) + ggtitle(variable) +
    scale_fill_manual(values=my_palette[1:5], name='') + 
    theme(legend.text = element_text(size=13),
          legend.title = element_text(size=15),
          axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          title=element_text(size=18))
  return(list(acp=acp_plot, var=var_plot, data=my_df))
}

####ABUNDANCE
fda_ab <- my_fda(data_ab, weight_ab, 10, 'Abondance')
#fda_ab$acp
#fda_ab$var

####NUTRIENTS
fda_nutri <- my_fda(data_nutri, weight_nutri, 5, 'Nutriments', nutri=TRUE)
#fda_nutri$acp
#fda_nutri$var

####DIFFUSION
fda_diff <- my_fda(data_diff, weight_diff, 5, 'Diffusion')
#fda_diff$acp
#fda_diff$var

