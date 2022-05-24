#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Functional PCA
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(ggpubr)
library(fda)

###source code###
source("src/my_palette.R")

###save preprocessed data###
save <- FALSE

###Preprocessing function for PCA###
build_ymat <- function(data, var){
  `%!in%` <- Negate(`%in%`)
  no_phyto <- c('HNABACC', 'LNABACC', 'TBACC', 'HNABACSSC', 'LNABACSSC', 'TBACSSC')
  data <- data %>% filter(GROUPE %!in% no_phyto)
  id <- which(colnames(data)==var)
  colnames(data)[id] <- 'VARIABLE'
  data$DATE <- ymd(data$DATE)
  data <- data %>% mutate("YEAR" = year(DATE), 
                          "MONTH"=month(DATE))
  moyennes <- data %>% 
    dplyr::group_by(YEAR, MONTH, GROUPE, ID_SITE) %>% 
    summarise("mean"=mean(log(VARIABLE+1)))
  moyennes <- na.omit(moyennes)
  moyennes <- as.data.frame(moyennes)
  years <- unique(moyennes$YEAR)
  gp <- unique(moyennes$GROUPE)
  sites <- unique(moyennes$ID_SITE)
  l <- length(sites)
  n <- length(years)
  m <- length(gp)
  ymat <- NULL
  curves <- NULL
  spp <- NULL
  curves_sites <- NULL
  for(s in 1:l){
    one_site <- moyennes %>% filter(ID_SITE==sites[s])
  for(i in 1:n){
    an <- one_site %>% filter(YEAR==years[i])
    for(j in 1:m){
      an_sp <- an %>% filter(GROUPE==gp[j])
      if(nrow(an_sp)==12){
        ymat <- cbind(ymat, an_sp$mean)
        curves <- c(curves, years[i])
        spp <- c(spp, gp[j])
        curves_sites <- c(curves_sites, sites[s])
    }}}
  }
  all_info <- data.frame("ANNEE"=curves, "SITE"=curves_sites, "GROUPE"=spp)
  return(list('ymat'=ymat, 'all_info'=all_info))
}

###FPCA, MEAN PERTURBATION PLOT, SCORE PLOT###
my_fpca <- function(y_mat, all_info, variable, sp=NULL, year_lab=FALSE, nutri=FALSE, norm=TRUE){
  sp_name <- data.frame(code=c("SYNC", "CRYC", "PROC", "PICOEC", "NANOEC"),
                        name=c("Synechococcus", "Cryptophytes", "Prochlorococcus",
                               "Pico-eucaryotes", "Nano-eucaryotes"))
  if(variable=='Diffusion'){
    sp_name <- data.frame(code=c("SYNSSC", "CRYSSC", "PROSSC", "PICOESSC", "NANOESSC"),
                          name=c("Synechococcus", "Cryptophytes", "Prochlorococcus",
                                 "Pico-eucaryotes", "Nano-eucaryotes"))
  }
  groupe <- 'Groupe'
  if(nutri==TRUE){
    sp_name <- data.frame(code=c("NH4", "NO3", "NO2", "PO4", "SIOH4"),
                                        name=c("Amonium", "Nitrate", "Nitrite",
                                               "Phosphate", "Silice"))
    groupe <- 'Nutriment'
  }
  y_mat <- as.matrix(y_mat)
  colnames(y_mat) <- NULL
  n_points <- nrow(y_mat)
  n_curves <- ncol(y_mat)
  arg_time <- seq(1, 12, 1)
  n_basis <- ceiling(arg_time/2)
  argvals <- matrix(arg_time, nrow=n_points, ncol=n_curves)
  #normalisation
  if(norm==TRUE){
  sum_mat <- apply(y_mat, 2, sum)
  y_mat <- sweep(y_mat, 2, sum_mat, FUN="/")
  }
  #create basis
  basis <- create.bspline.basis(c(1, 12), nbasis=n_basis, norder=4)
  l_use <- 1e-4 #lambda
  #select sp
  if(is.null(sp)==FALSE){
    ind <- which(all_info$GROUPE %in% sp)
    argvals <- argvals[,ind]
    y_mat <- y_mat[,ind]
    all_info <- all_info[ind,]
    sp_name <- sp_name %>% filter(code %in% sp)
  }
  W.obj <- Data2fd(argvals = argvals, y = y_mat, basisobj = basis, lambda = l_use)
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
  xx <- seq(1, 12, 0.2)
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
  pc3 <- matrix(facp$harmonics$coefs[,3], ncol=1)
  pert_3p$coefs <- mean$coefs + sqrt(facp$values[3]) * pc3
  pert_3n$coefs <- mean$coefs - sqrt(facp$values[3]) * pc3
  df_pcp3 <- data.frame(x=xx, yhat = predict(pert_3p, newdata=xx))
  df_pcn3 <- data.frame(x=xx, yhat = predict(pert_3n, newdata=xx))
  mois <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
            "Oct", "Nov", "Dec")
  A <- ggplot() + geom_path(data=df_mean, aes(xx, mean), cex=1) +
    geom_point(data=df_pcp1, aes(xx, mean), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn1, aes(xx, mean), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC1', pc1_pi, '%'))
  B <- ggplot() + geom_path(data=df_mean, aes(xx, mean), cex=1) +
    geom_point(data=df_pcp2, aes(xx, mean), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn2, aes(xx, mean), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC2', pc2_pi, '%'))
  C <- ggplot() + geom_path(data=df_mean, aes(xx, mean), cex=1) +
    geom_point(data=df_pcp3, aes(xx, mean), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn3, aes(xx, mean), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC3', pc3_pi, '%'))
  plot_pert <- ggarrange(A,B, C, ncol=1, nrow=3)
  ##Components scores
  df_pc <- data.frame('PC1'=facp$scores[,1], 
                      'PC2'=facp$scores[,2], 
                      'PC3'=facp$scores[,3])
  my_df <- cbind(df_pc, all_info)
  my_df$SITE <- factor(my_df$SITE, 
                       labels=c('Banyuls', 'Marseille', 'Villefranche'))
  my_df$GROUPE <- factor(my_df$GROUPE, 
                         levels=sp_name$code,
                         labels=sp_name$name)
  plot_score <- my_df %>% ggplot() + 
    theme_light() + 
    geom_vline(xintercept = 0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed')  +
    theme(aspect.ratio = 1) +
    xlab(paste('PC1', pc1_pi, '%')) +
    ylab(paste('PC2', pc2_pi, '%')) +
    scale_color_manual(values=my_palette[2:4], name='Site') +
    theme(legend.text = element_text(size=13),
          legend.title = element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  if(length(sp)==1){
    solo <- unique(my_df$GROUPE)
    plot_score <- plot_score + 
      geom_point(aes(PC1, PC2, col=SITE), size=2) +
      ggtitle(solo)
  }
  if(length(sp)!=1){
    plot_score <- plot_score +
      geom_point(aes(PC1, PC2, col=SITE, shape=GROUPE), size=2) +
      scale_shape_discrete(name=groupe) 
  }
  if(year_lab==TRUE){
    plot_score <- plot_score +
      geom_text(aes(PC1, PC2, label=ANNEE, col=factor(SITE)), 
                vjust="inward", hjust="inward")}
  my_plot <- ggarrange(plot_pert, NULL, plot_score, nrow=1, 
                       widths=c(1, 0.1, 1.2))
  return(my_plot)
}  

###load raw data###
data_ab <- read.csv("data/PICONANO_AB.csv")
data_nutri <- read.csv("data/NUTRIENTS.csv")
data_diff <- read.csv("data/PICONANO_DIFF.csv")

###save preprocessed data###
if(save==TRUE){
  export <- build_ymat(data_ab, 'ABONDANCE')
  write.csv(export$ymat, 'results/FPCA_YMAT_AB.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/FPCA_INFO_AB.csv', row.names=FALSE)
  export <- build_ymat(data_nutri, 'CONCENTRATION')
  write.csv(export$ymat, 'results/FPCA_YMAT_NUTRI.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/FPCA_INFO_NUTRI.csv', row.names=FALSE)
  export <- build_ymat(data_diff, 'DIFFUSION')
  write.csv(export$ymat, 'results/FPCA_YMAT_DIFF.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/FPCA_INFO_DIFF.csv', row.names=FALSE)
}

###load preprocess data###
ymat_ab <- read.csv('results/FPCA_YMAT_AB.csv')
info_ab <- read.csv('results/FPCA_INFO_AB.csv')
ymat_nutri <- read.csv('results/FPCA_YMAT_NUTRI.csv')
info_nutri <- read.csv('results/FPCA_INFO_NUTRI.csv')
ymat_diff <- read.csv('results/FPCA_YMAT_DIFF.csv')
info_diff <- read.csv('results/FPCA_INFO_DIFF.csv')

###fpca on abundance###
out1 <- which(info_ab$GROUPE=='PROC' & 
                    info_ab$SITE==12 & 
                    info_ab$ANNEE==2016)
out2 <- which(info_ab$GROUPE=='PICOEC' & 
                info_ab$SITE==11 & 
                info_ab$ANNEE==2013)
outs <- c(out1, out2)
ymat_ab <- ymat_ab[,-outs]
info_ab <- info_ab[-outs,]
fpca_ab <- my_fpca(ymat_ab, info_ab, 'Abondance')
fpca_cry <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='CRYC', year_lab=TRUE, norm=FALSE)
fpca_pro <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='PROC', year_lab=TRUE, norm=FALSE)
fpca_syn <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='SYNC', year_lab=TRUE, norm=FALSE)
fpca_picoe <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='PICOEC', year_lab=TRUE, norm=FALSE)
fpca_nanoe <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='NANOEC', year_lab=TRUE, norm=FALSE)
fpca_other <- my_fpca(ymat_ab, info_ab, 'Abondance', sp=c('SYNC', 'PICOEC', 'NANOEC'))
fpca_picoe
fpca_cry
fpca_other


###fpca on nutrients###
out1 <- which(info_nutri$GROUPE=='NH4' & 
                info_nutri$SITE==11 & 
                info_nutri$ANNEE==2018)
out2 <- which(info_nutri$GROUPE=='PO4' & 
                info_nutri$SITE==12 & 
                info_nutri$ANNEE==2014)
out3 <- which(info_nutri$GROUPE=='PO4' & 
                info_nutri$SITE==11 & 
                info_nutri$ANNEE==2014)
outs <- c(out1, out2, out3)
ymat_nutri <- ymat_nutri[,-outs]
info_nutri <- info_nutri[-outs,]
fpca_nutri <- my_fpca(ymat_nutri, info_nutri, 'Concentration', nutri=TRUE)
fpca_N <- my_fpca(ymat_nutri, info_nutri, 'Concentration', 
                      sp=c('NH4', 'NO3', 'NO2'), 
                      nutri=TRUE)
fpca_NH4 <- my_fpca(ymat_nutri, info_nutri, 'Concentration', sp='NH4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NO3 <- my_fpca(ymat_nutri, info_nutri, 'Concentration', sp='NO3', 
                    year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NO2 <- my_fpca(ymat_nutri, info_nutri, 'Concentration', sp='NO2', 
                    year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_P <- my_fpca(ymat_nutri, info_nutri, 'Concentration', sp='PO4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_S <- my_fpca(ymat_nutri, info_nutri, 'Concentration', sp='SIOH4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NH4
fpca_NO3
fpca_P

###fpca on diffusion###
out1 <- which(info_diff$GROUPE=='PICOESSC' & 
                info_diff$SITE==11 & 
                info_diff$ANNEE==2019)
out2 <- which(info_diff$GROUPE=='NANOESSC' & 
                info_diff$SITE==11 & 
                info_diff$ANNEE==2013)
outs <- c(out1, out2)
ymat_diff <- ymat_diff[,-outs]
info_diff <- info_diff[-outs,]
fpca_diff <- my_fpca(ymat_diff, info_diff, 'Diffusion')
fpca_dcry <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='CRYSSC', year_lab=TRUE, norm=FALSE)
fpca_dsyn <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='SYNSSC', year_lab=TRUE, norm=FALSE)
fpca_dpicoe <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='PICOESSC', year_lab=TRUE, norm=FALSE)
fpca_dnanoe <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='NANOESSC', year_lab=TRUE, norm=FALSE)

