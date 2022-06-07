#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Functional PCA CTD
################

###working directory###
setwd("C:/Users/ambre/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(ggpubr)
library(fda)

###source code###
source("src/my_palette.R")

###preprocessing data and lambda choice function###
pre_fpca <- function(data, pmin=0.5, pmax=55, l_use){
  #interpolation for some prof
  arg_prof <- seq(pmin, pmax, 0.5)
  n_points <- length(arg_prof)
  dates <- unique(data$DATE)
  all_temp <- data.frame("DATE"=NULL, "arg_prof"=NULL, "temp"=NULL)
  for(d in dates){
    profil <- data %>% filter(DATE==d)
    temp <- approx(profil$PROFONDEUR, profil$TEMPERATURE, rule=2, xout=arg_prof)$y
    date <- rep(d, n_points)
    seul <- data.frame("DATE"=date, arg_prof, temp)
    all_temp <- rbind(all_temp, seul)
  }
  all_temp <- all_temp %>% mutate("DATE"=ymd(DATE), 
                                  "ANNEE"= year(DATE), "MOIS"=month(DATE), 
                                  "JOUR"=mday(DATE))
  #cubic spline smoothing
  n_basis <- floor(n_points/2)
  n_curves <- length(dates)
  argvals <- matrix(all_temp$arg_prof, nrow=n_points, ncol=n_curves)
  y_mat <- matrix(all_temp$temp, nrow=n_points, ncol=n_curves)
  basis <- create.bspline.basis(c(0.5,55), nbasis=n_basis, norder=4)
  ##ajust lambda by vis 8 curves
  set.seed(2)
  ind <- sample.int(n_curves, 8)
  par(mar=c(2.1,2.3,1.5,0.5), mgp=c(1.2,0.5,0))
  W.obj_test <- Data2fd(argvals = argvals[,ind], y = y_mat[,ind], basisobj = basis, lambda = l_use)
  xx <- seq(pmin, pmax, 0.2)
  df_all <- NULL
  for(i in 1:8){
    prof <- fd(basisobj=basis)
    prof$coefs <- W.obj_test$coefs[,i]
    df_prof <- data.frame(N=i, x=xx, yhat = predict(prof, newdata=xx))
    df_all <- rbind(df_all, df_prof)
  }
  original <- data.frame(x=c(argvals[,ind]), yhat=c(y_mat[,ind]))
  l_test <- df_all %>% ggplot()  +
    geom_path(aes(yhat, x, group=N), size=.7) +
    geom_point(data=original, aes(yhat, x), shape=1, col=my_palette[2], size=3) +
    scale_y_reverse(breaks=seq(0, 50, 10)) +
    theme_light() + ylab('Profondeur (m)') + xlab('Température (°C)') +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18))
  W.obj <- Data2fd(argvals = argvals, y = y_mat, basisobj = basis, lambda = l_use)
  label <- unique(all_temp$DATE)
  return(list(plot= l_test, W.obj=W.obj, dates=label, basis=basis, df_test=df_all, pt_test=original))
}

###FPCA, MEAN PERTURBATION PLOT, SCORE PLOT###
my_fpca <- function(W.obj, basis, label, id_site, pmin=0.5, pmax=50){
  ##FACP
  facp <- pca.fd(fdobj = W.obj, nharm=5)
  ##Components indertia
  pc1_pi <- round(facp$values[1]/sum(facp$values)*100, 2)
  pc2_pi <- round(facp$values[2]/sum(facp$values)*100, 2)
  pc3_pi <- round(facp$values[3]/sum(facp$values)*100, 2)
  ##PC as mean perturbation
  pert_1p <- fd(basisobj=basis)
  pert_1n <- fd(basisobj=basis)
  pert_2p <- fd(basisobj=basis)
  pert_2n <- fd(basisobj=basis)
  pert_3p <- fd(basisobj=basis)
  pert_3n <- fd(basisobj=basis)
  xx <- seq(pmin, pmax, 0.7)
  mean <- mean.fd(W.obj)
  pc1 <- matrix(facp$harmonics$coefs[,1], ncol=1)
  pert_1p$coefs <- mean$coefs + sqrt(facp$values[1]) * pc1
  pert_1n$coefs <- mean$coefs - sqrt(facp$values[1]) * pc1
  df_mean <- data.frame(x=xx, yhat = predict(mean, newdata=xx))
  df_pcp1 <- data.frame(x=xx, yhat = predict(pert_1p, newdata=xx))
  df_pcn1 <- data.frame(x=xx, yhat = predict(pert_1n, newdata=xx))
  pc2 <- matrix(facp$harmonics$coefs[,2], ncol=1)
  pert_2p$coefs <- mean$coefs + sqrt(facp$values[2]) * pc2
  pert_2n$coefs <- mean$coefs - sqrt(facp$values[2]) * pc2
  df_pcp2 <- data.frame(x=xx, yhat = predict(pert_2p, newdata=xx))
  df_pcn2 <- data.frame(x=xx, yhat = predict(pert_2n, newdata=xx))
  A <- ggplot() + geom_path(data=df_mean, aes(mean, xx), cex=1) +
    geom_point(data=df_pcp1, aes(mean, xx), cex=1.5, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn1, aes(mean, xx), cex=1, col=my_palette[2], linetype='dashed') +
    scale_y_reverse() +
    theme_minimal() + ylab('Profondeur (m)') + xlab('Température (°C)') +
    ggtitle(paste('PC1', pc1_pi, '%')) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          title = element_text(size=20))
  B <- ggplot() + geom_path(data=df_mean, aes(mean, xx), cex=1) +
    geom_point(data=df_pcp2, aes(mean, xx), cex=1.5, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn2, aes(mean, xx), cex=1, col=my_palette[2], linetype='dashed') +
    scale_y_reverse() +
    theme_minimal() + ylab('') + xlab('Température (°C)') +
    ggtitle(paste('PC2', pc2_pi, '%')) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          title = element_text(size=20))
  plot_pert <- ggarrange(A, NULL, B, ncol=3, nrow=1, widths=c(1, .2, 1))
  df_pc <- data.frame('DATE'=label,
                      'ID_SITE'=id_site,
    'PC1'=facp$scores[,1], 
    'PC2'=facp$scores[,2])
  return(list(plot=plot_pert, scores=df_pc))
}  


###load data###
ctd_mar <- read.csv("data/CTD_MARSEILLE.csv")

###FPCA marseille###
dat_mar <- pre_fpca(ctd_mar, l_use=1e-1)
fpca_mar <- my_fpca(dat_mar$W.obj, dat_mar$basis, dat_mar$dates, 11)
M <- ggarrange(NULL, fpca_mar$plot, NULL, ncol=3, nrow=1, widths=c(.1, 1, .1))
export <- fpca_mar$scores
#write.csv(export, "results/DATA_PC_CTD_MARSEILLE.csv", row.names=FALSE)

export <- export %>% mutate('col1'=ifelse(PC1>0, 1, 2), 'col2'=ifelse(PC2>0, 1, 2))
date_range <- ymd(c("1995-01-01", "2022-01-01"))
A <- export %>% ggplot() + 
  geom_segment(aes(x=DATE, y=0, xend=DATE, yend=PC1, col=factor(col1)), size=.7) +
  scale_x_date(limits=date_range, 
               breaks=seq.Date(ymd("1995-01-01"), ymd("2022-01-01"), by='2 year'),
               labels=seq(1995, 2022, 2),
               expand=c(1e-2,1e-2)) +
  scale_color_manual(values=my_palette, name='') +
  theme_light() + ylab('PC1') + xlab('') + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.position='none')
B <- export %>% ggplot() + 
  geom_segment(aes(x=DATE, y=0, xend=DATE, yend=PC2, col=factor(col2)), size=.7) +
  scale_x_date(limits=date_range, 
               breaks=seq.Date(ymd("1995-01-01"), ymd("2022-01-01"), by='2 year'),
               labels=seq(1995, 2022, 2),
               expand=c(1e-2,1e-2)) +
  scale_color_manual(values=my_palette, name='') +
  theme_light() + ylab('PC2') + xlab('Temps') + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18),
        legend.position='none')

S <- ggarrange(A, B, ncol=1, nrow=2)
ggarrange(M, NULL, S, ncol=1, nrow=3, heights=c(1, .1, 1.1), labels=c('a)', '', 'b)'))

