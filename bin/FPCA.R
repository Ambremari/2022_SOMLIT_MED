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
build_ymat <- function(data, site, sp){
  data$DATE <- ymd(data$DATE)
  colnames(data)[3] <- 'var'
  data <- data %>% mutate("YEAR" = year(DATE), 
                          "MONTH"=month(DATE), 
                          "HALF"=ifelse(mday(DATE)<=15, 1, 2))
  moyennes <- data %>% 
    dplyr::group_by(YEAR, MONTH, HALF) %>% 
    summarise("mean"=mean(var))
  moyennes <- as.data.frame(moyennes)
  years <- unique(moyennes$YEAR)
  n <- length(years)
  ymat <- NULL
  curves <- NULL
  for(i in 1:n){
    an <- moyennes %>% filter(YEAR==years[i])
    first_line <- head(an, 1)
    if(first_line$MONTH==1 & first_line$HALF==2){
      debut <- first_line %>% mutate(HALF=1)
      an <- rbind(debut, an)}
    last_line <- tail(an, 1)
    if(last_line$MONTH==12 & last_line$HALF==1){
      fin <- last_line %>% mutate(HALF=2)
      an <- rbind(an, fin)}
    if(nrow(an)==24){
      ymat <- cbind(ymat, an$mean)
      curves <- c(curves, years[i])
      n_curves <- length(curves)
      sites <- rep(site, n_curves)
      spp <- rep(sp, n_curves)}
  }
  return(list('ymat'=ymat, 'years'=curves, 'site'=sites, 'sp'=spp))
}

merge_mat <- function(data, site, sp, y_mat, all_info, nutri=FALSE){
  single_ts <- build_ymat(data, site, sp)
  new_mat <- cbind(y_mat, single_ts$ymat)
  single_info <- data.frame("ANNEE"=single_ts$years, "SITE"=single_ts$site, "GROUPE"=single_ts$sp)
  new_info <- rbind(all_info, single_info)
  return(list('y_mat'=new_mat, 'all_info'=new_info))
}

###FPCA, MEAN PERTURBATION PLOT, SCORE PLOT###
my_fpca <- function(y_mat, all_info, variable, sp=NULL, year_lab=FALSE, nutri=FALSE, norm=TRUE){
  sp_name <- data.frame(code=c("SYN", "CRY", "PRO", "PICOE", "NANOE"),
                        name=c("Synechococcus", "Cryptophytes", "Prochlorococcus",
                               "Pico-eucaryotes", "Nano-eucaryotes"))
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
  arg_time <- seq(1, 12.5, 0.5)
  argvals <- matrix(arg_time, nrow=n_points, ncol=n_curves)
  #normalisation
  if(norm==TRUE){
  sum_mat <- apply(y_mat, 2, sum)
  y_mat <- sweep(y_mat, 2, sum_mat, FUN="/")
  }
  #create basis
  basis <- create.bspline.basis(c(1,12.5), breaks=arg_time, norder=4)
  l_use <- 0.005 #lambda
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
  xx <- seq(1, 12.5, 0.2)
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
    geom_path(data=df_pcp1, aes(xx, mean), cex=1, col=my_palette[1], linetype='dashed') +
    geom_path(data=df_pcn1, aes(xx, mean), cex=1, col=my_palette[2], linetype='dotted') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC1', pc1_pi, '%'))
  B <- ggplot() + geom_path(data=df_mean, aes(xx, mean), cex=1) +
    geom_path(data=df_pcp2, aes(xx, mean), cex=1, col=my_palette[1], linetype='dashed') +
    geom_path(data=df_pcn2, aes(xx, mean), cex=1, col=my_palette[2], linetype='dotted') +
    theme_light() + ylab(variable) + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC2', pc2_pi, '%'))
  C <- ggplot() + geom_path(data=df_mean, aes(xx, mean), cex=1) +
    geom_path(data=df_pcp3, aes(xx, mean), cex=1, col=my_palette[1], linetype='dashed') +
    geom_path(data=df_pcn3, aes(xx, mean), cex=1, col=my_palette[2], linetype='dotted') +
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


###MERGE ALL ABUNDANCE DATA###
if(save==TRUE){
  y_mat <- NULL
  all_info <- NULL
  data <- read.csv("results/TS_SYN_AB_BANYULS.csv")
  merged <- merge_mat(data, 10, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SYN_AB_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SYN_AB_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_AB_BANYULS.csv")
  merged <- merge_mat(data, 10, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_AB_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_AB_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_AB_BANYULS.csv")
  merged <- merge_mat(data, 10, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_AB_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_AB_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_AB_BANYULS.csv")
  merged <- merge_mat(data, 10, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_AB_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_AB_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_AB_BANYULS.csv")
  merged <- merge_mat(data, 10, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_AB_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_AB_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  write.csv(y_mat, "results/FPCA_YMAT_AB.csv", row.names=FALSE)
  write.csv(all_info, "results/FPCA_INFO_AB.csv", row.names=FALSE)
}

###MERGE ALL DIFFUSION DATA###
if(save==TRUE){
  y_mat <- NULL
  all_info <- NULL
  data <- read.csv("results/TS_SYN_DIFF_BANYULS.csv")
  merged <- merge_mat(data, 10, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SYN_DIFF_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SYN_DIFF_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "SYN", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_DIFF_BANYULS.csv")
  merged <- merge_mat(data, 10, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_DIFF_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_CRY_DIFF_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "CRY", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_DIFF_BANYULS.csv")
  merged <- merge_mat(data, 10, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_DIFF_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PRO_DIFF_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "PRO", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_DIFF_BANYULS.csv")
  merged <- merge_mat(data, 10, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_DIFF_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PICOE_DIFF_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "PICOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_DIFF_BANYULS.csv")
  merged <- merge_mat(data, 10, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_DIFF_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NANOE_DIFF_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "NANOE", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  write.csv(y_mat, "results/FPCA_YMAT_DIFF.csv", row.names=FALSE)
  write.csv(all_info, "results/FPCA_INFO_DIFF.csv", row.names=FALSE)
}


###MERGE ALL NUTRIENTS DATA###
if(save==TRUE){
  y_mat <- NULL
  all_info <- NULL
  data <- read.csv("results/TS_NH4_BANYULS.csv")
  merged <- merge_mat(data, 10, "NH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NH4_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "NH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NH4_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "NH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO3_BANYULS.csv")
  merged <- merge_mat(data, 10, "NO3", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO3_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "NO3", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO3_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "NO3", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO2_BANYULS.csv")
  merged <- merge_mat(data, 10, "NO2", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO2_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "NO2", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_NO2_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "NO2", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PO4_BANYULS.csv")
  merged <- merge_mat(data, 10, "PO4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PO4_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "PO4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_PO4_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "PO4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SIOH4_BANYULS.csv")
  merged <- merge_mat(data, 10, "SIOH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SIOH4_MARSEILLE.csv")
  merged <- merge_mat(data, 11, "SIOH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  data <- read.csv("results/TS_SIOH4_VILLEFRANCHE.csv")
  merged <- merge_mat(data, 12, "SIOH4", y_mat, all_info)
  y_mat <- merged$y_mat
  all_info <- merged$all_info
  write.csv(y_mat, "results/FPCA_YMAT_NUTRI.csv", row.names=FALSE)
  write.csv(all_info, "results/FPCA_INFO_NUTRI.csv", row.names=FALSE)
}

###load preprocessed data###
y_mat_ab <- read.csv("results/FPCA_YMAT_AB.csv")
all_info_ab <- read.csv("results/FPCA_INFO_AB.csv")
y_mat_diff <- read.csv("results/FPCA_YMAT_DIFF.csv")
all_info_diff <- read.csv("results/FPCA_INFO_DIFF.csv")
y_mat_nutri <- read.csv("results/FPCA_YMAT_NUTRI.csv")
all_info_nutri <- read.csv("results/FPCA_INFO_NUTRI.csv")

###fpca on abundance###
fpca_ab <- my_fpca(y_mat_ab, all_info_ab, 'Abondance')
fpca_cry <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp='CRY', year_lab=TRUE)
fpca_pro <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp='PRO', year_lab=TRUE)
fpca_syn <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp='SYN', year_lab=TRUE, norm=FALSE)
fpca_picoe <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp='PICOE', year_lab=TRUE)
fpca_nanoe <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp='NANOE', year_lab=TRUE)
fpca_other <- my_fpca(y_mat_ab, all_info_ab, 'Abondance', sp=c('SYN', 'PICOE', 'NANOE'))

###fpca on nutrients###
fpca_nutri <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', nutri=TRUE)
fpca_N <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', 
                      sp=c('NH4', 'NO3', 'NO2'), 
                      nutri=TRUE)
fpca_NH4 <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', sp='NH4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NO3 <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', sp='NO3', 
                    year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NO2 <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', sp='NO2', 
                    year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_P <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', sp='PO4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_S <- my_fpca(y_mat_nutri, all_info_nutri, 'Concentration', sp='SIOH4', 
                  year_lab=TRUE, nutri=TRUE, norm=FALSE)
fpca_NH4
fpca_NO3
fpca_P

###fpca on diffusion###
pb_pro <- which(all_info_diff$SITE==10 &
                all_info_diff$GROUPE=='PRO' & 
                all_info_diff$ANNEE<2015)
y_mat_diff <- y_mat_diff[,-pb_pro]
all_info_diff <- all_info_diff[-pb_pro,]
fpca_diff <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion')
fpca_dcry <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp='CRY', year_lab=TRUE, norm=FALSE)
fpca_dpro <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp='PRO', year_lab=TRUE, norm=FALSE)
fpca_dsyn <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp='SYN', year_lab=TRUE, norm=FALSE)
fpca_dpicoe <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp='PICOE', year_lab=TRUE, norm=FALSE)
fpca_dnanoe <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp='NANOE', year_lab=TRUE, norm=FALSE)
fpca_dother <- my_fpca(y_mat_diff, all_info_diff, 'Diffusion', sp=c('SYN', 'PRO'), norm=FALSE)

