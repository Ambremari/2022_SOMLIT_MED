#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Functional PCA par
################

###working directory###
setwd("C:/Users/ambre/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(ggpubr)
library(fdapace)

###source code###
source("src/my_palette.R")


###Preprocessing function for PCA###
build_ymat <- function(data, var){
  id <- which(colnames(data)==var)
  colnames(data)[id] <- 'VARIABLE'
  data$DATE <- ymd(data$DATE)
  data <- data %>% mutate("YEAR" = year(DATE), 
                          "MONTH"=month(DATE))
  moyennes <- data %>% 
    dplyr::group_by(YEAR, MONTH, ID_SITE) %>% 
    summarise("mean"=mean(VARIABLE))
  moyennes <- na.omit(moyennes)
  moyennes <- as.data.frame(moyennes)
  years <- unique(moyennes$YEAR)
  sites <- unique(moyennes$ID_SITE)
  l <- length(sites)
  n <- length(years)
  ymat <- NULL
  curves <- NULL
  curves_sites <- NULL
  for(s in 1:l){
    one_site <- moyennes %>% filter(ID_SITE==sites[s])
    for(i in 1:n){
      an <- one_site %>% filter(YEAR==years[i])
        my_month <- unique(an$MONTH)
        my_means <- NULL
        for(z in 1:12){
          my_mean <- ifelse(z %in% my_month, an$mean[z], NA)
          my_means <- rbind(my_means, my_mean)
        }
        if(sum(is.na(my_means))!=12){
          ymat <- cbind(ymat, my_means)
          curves <- c(curves, years[i])
          curves_sites <- c(curves_sites, sites[s])
        }
    }}
  all_info <- data.frame("ANNEE"=curves, "SITE"=curves_sites)
  return(list('ymat'=ymat, 'all_info'=all_info))
}

###FPCA, MEAN PERTURBATION PLOT, SCORE PLOT###
my_fpca <- function(y_mat, all_info, variable, year_lab=FALSE){
  y_mat <- as.matrix(y_mat)
  colnames(y_mat) <- NULL
  n_points <- nrow(y_mat)
  n_curves <- ncol(y_mat)
  arg_time <- seq(1, 12, 1)
  argvals <- matrix(arg_time, nrow=n_points, ncol=n_curves)
  ##FACP
  facp <- FPCA(data.frame(y_mat), data.frame(argvals),
               optns=list(useBinnedData='OFF',
                          methodMuCovEst='smooth',
                          userBwCov=3,
                          userBwMu=2,
                          kernel="quar",
                          nRegGrid=100,
                          dataType="Sparse",
                          methodXi='CE'
               ))
  
  ##Components inertia
  pc1_pi <- round(facp$lambda[1]/sum(facp$lambda)*100, 2)
  pc2_pi <- round(facp$lambda[2]/sum(facp$lambda)*100, 2)
  ##Components weights for fda
  k2 <- sqrt(facp$lambda[2]/facp$lambda[1])
  weights <- c(1, k2)
  ##PC as mean perturbation
  xx <- facp$workGrid
  mean <- facp$mu
  pc1 <- facp$phi[,1]
  pert_1p <- mean + sqrt(facp$lambda[1]) * pc1
  pert_1n <- mean - sqrt(facp$lambda[1]) * pc1
  df_mean <- data.frame(x=xx, yhat = mean)
  df_pcp1 <- data.frame(x=xx, yhat = pert_1p)
  df_pcn1 <- data.frame(x=xx, yhat = pert_1n)
  pc2 <- facp$phi[,2]
  pert_2p <- mean + sqrt(facp$lambda[2]) * pc2
  pert_2n <- mean - sqrt(facp$lambda[2]) * pc2
  df_pcp2 <- data.frame(x=xx, yhat = pert_2p)
  df_pcn2 <- data.frame(x=xx, yhat = pert_2n)
  mois <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
            "Oct", "Nov", "Dec")
  A <- ggplot() + geom_path(data=df_mean, aes(xx, yhat), cex=1) +
    geom_point(data=df_pcp1, aes(xx, yhat), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn1, aes(xx, yhat), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC1', pc1_pi, '%'))+
    theme( axis.text=element_text(size=16),
           axis.title=element_text(size=18),
           title=element_text(size=18))
  B <- ggplot() + geom_path(data=df_mean, aes(xx, yhat), cex=1) +
    geom_point(data=df_pcp2, aes(xx, yhat), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn2, aes(xx, yhat), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC2', pc2_pi, '%'))+
    theme( axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          title=element_text(size=18))
  plot_pert <- ggarrange(A, NULL, B, ncol=1, nrow=3, heights=c(1, .1, 1))
  ##Components scores
  df_pc <- data.frame('PC1'=facp$xiEst[,1], 
                      'PC2'=facp$xiEst[,2])
  my_df <- cbind(df_pc, all_info)
  my_df$SITE <- factor(my_df$SITE,
                       levels=c(10, 11, 12),
                       labels=c('Banyuls', 'Marseille', 'Villefranche'))
  plot_score <- my_df %>% ggplot() + 
    theme_light() + 
    geom_vline(xintercept = 0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed')  +
    geom_point(aes(PC1, PC2, col=ANNEE), size=2) +
    scale_colour_gradient(low=my_palette[2], high=my_palette[3], name='Date') +
    xlab(paste('PC1', pc1_pi, '%')) +
    ylab(paste('PC2', pc2_pi, '%')) +
    theme(aspect.ratio = 1,
          legend.text = element_text(size=16),
          legend.title = element_text(size=18),
          axis.text=element_text(size=16),
          axis.title=element_text(size=18))
  if(year_lab==TRUE){
    plot_score <- plot_score +
      geom_text(aes(PC1, PC2, label=ANNEE, col=ANNEE), 
                vjust="inward", hjust="inward", size=7)}
  my_plot <- ggarrange(plot_pert, NULL, plot_score, nrow=1, 
                       widths=c(1, 0.1, 1.2))
  return(list("data"=my_df, "plot"=my_plot, 'weights'=weights, 'PC1'=A, 'PC2'=B))
} 

###load raw data###
data_ctd <- read.csv("results/DATA_PC_CTD_MARSEILLE.csv")


###PC1###
mat_pc1 <- build_ymat(data_ctd, 'PC1')
fpca_pc1 <- my_fpca(mat_pc1$ymat, mat_pc1$all_info, 'PC1', year_lab=TRUE)
fpca_pc1$plot

###PC2
mat_pc2 <- build_ymat(data_ctd, 'PC2')
fpca_pc2 <- my_fpca(mat_pc2$ymat, mat_pc2$all_info, 'PC2', year_lab=TRUE)
fpca_pc2$plot

