#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Functional PCA par
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(ggpubr)
library(fdapace)

###source code###
source("src/my_palette.R")

###save preprocess data
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
        my_month <- unique(an_sp$MONTH)
        my_means <- NULL
        for(z in 1:12){
          my_mean <- ifelse(z %in% my_month, an_sp$mean[z], NA)
          my_means <- rbind(my_means, my_mean)
        }
        if(sum(is.na(my_means))!=12){
          ymat <- cbind(ymat, my_means)
          curves <- c(curves, years[i])
          spp <- c(spp, gp[j])
          curves_sites <- c(curves_sites, sites[s])
        }}
    }}
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
  argvals <- matrix(arg_time, nrow=n_points, ncol=n_curves)
  #normalisation
  if(norm==TRUE){
    sum_mat <- apply(y_mat, 2, sum, na.rm=TRUE)
    y_mat <- sweep(y_mat, 2, sum_mat, FUN="/")
  }
  #select sp
  if(is.null(sp)==FALSE){
    ind <- which(all_info$GROUPE %in% sp)
    argvals <- argvals[,ind]
    y_mat <- y_mat[,ind]
    all_info <- all_info[ind,]
    sp_name <- sp_name %>% filter(code %in% sp)
  }
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
  pc3_pi <- round(facp$lambda[3]/sum(facp$lambda)*100, 2)
  inertia <- c(pc1_pi, pc2_pi, pc3_pi)
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
  pc3 <- facp$phi[,3]
  pert_3p <- mean + sqrt(facp$lambda[3]) * pc3
  pert_3n <- mean - sqrt(facp$lambda[3]) * pc3
  df_pcp3 <- data.frame(x=xx, yhat = pert_3p)
  df_pcn3 <- data.frame(x=xx, yhat = pert_3n)
  mois <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
            "Oct", "Nov", "Dec")
  A <- ggplot() + geom_path(data=df_mean, aes(xx, yhat), cex=1) +
    geom_point(data=df_pcp1, aes(xx, yhat), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn1, aes(xx, yhat), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab('variable') + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC1', pc1_pi, '%'))
  B <- ggplot() + geom_path(data=df_mean, aes(xx, yhat), cex=1) +
    geom_point(data=df_pcp2, aes(xx, yhat), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn2, aes(xx, yhat), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab('variable') + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC2', pc2_pi, '%'))
  C <- ggplot() + geom_path(data=df_mean, aes(xx, yhat), cex=1) +
    geom_point(data=df_pcp3, aes(xx, yhat), cex=1, col=my_palette[1], shape=3) +
    geom_path(data=df_pcn3, aes(xx, yhat), cex=1, col=my_palette[2], linetype='dashed') +
    theme_light() + ylab('variable') + xlab('Temps') +
    scale_x_continuous(breaks = 1:12, labels=mois) +
    ggtitle(paste('PC3', pc3_pi, '%'))
  plot_pert <- ggarrange(A,B, C, ncol=1, nrow=3)
  ##Components scores
  df_pc <- data.frame('PC1'=facp$xiEst[,1], 
                      'PC2'=facp$xiEst[,2], 
                      'PC3'=facp$xiEst[,3])
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
      geom_point(aes(PC1, PC2, col=factor(SITE), shape=GROUPE), size=2) 
      scale_shape_discrete(name=groupe) 
  }
  if(year_lab==TRUE){
    plot_score <- plot_score +
      geom_text(aes(PC1, PC2, label=ANNEE, col=factor(SITE)), 
                vjust="inward", hjust="inward")}
  my_plot <- ggarrange(plot_pert, NULL, plot_score, nrow=1, 
                       widths=c(1, 0.1, 1.2))
  return(list("data"=my_df, "plot"=my_plot, "inertia"=inertia))
}  

###load raw data###
data_ab <- read.csv("data/PICONANO_AB.csv")
data_nutri <- read.csv("data/NUTRIENTS.csv")
data_diff <- read.csv("data/PICONANO_DIFF.csv")

###save preprocessed data###
if(save==TRUE){
  export <- build_ymat(data_ab, 'ABONDANCE')
  write.csv(export$ymat, 'results/fpca_ymat_ab.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/fpca_info_ab.csv', row.names=FALSE)
  export <- build_ymat(data_nutri, 'CONCENTRATION')
  write.csv(export$ymat, 'results/fpca_ymat_nutri.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/fpca_info_nutri.csv', row.names=FALSE)
  export <- build_ymat(data_diff, 'DIFFUSION')
  write.csv(export$ymat, 'results/fpca_ymat_diff.csv', row.names=FALSE)
  write.csv(export$all_info, 'results/fpca_info_diff.csv', row.names=FALSE)
}

###load preprocess data###
ymat_ab <- read.csv('results/fpca_ymat_ab.csv')
info_ab <- read.csv('results/fpca_info_ab.csv')
ymat_nutri <- read.csv('results/fpca_ymat_nutri.csv')
info_nutri <- read.csv('results/fpca_info_nutri.csv')
ymat_diff <- read.csv('results/fpca_ymat_diff.csv')
info_diff <- read.csv('results/fpca_info_diff.csv')

###fpca on abundance###
out1 <- which(info_ab$ANNEE==2021)
out2 <- which(info_ab$SITE==12 & 
                info_ab$ANNEE==2019)
out3 <- which(info_ab$GROUPE=='PROC' & 
                info_ab$SITE==12 & 
                info_ab$ANNEE==2016)
out4 <- which(info_ab$GROUPE=='PICOEC' & 
                info_ab$SITE==11 & 
                info_ab$ANNEE==2013)
outs <- c(out1, out2, out3, out4)
ymat_ab <- ymat_ab[,-outs]
info_ab <- info_ab[-outs,]
fpca_ab <- my_fpca(ymat_ab, info_ab, 'Abondance')
fpca_cry <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='CRYC', year_lab=TRUE, norm=FALSE)
fpca_pro <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='PROC', year_lab=TRUE, norm=FALSE)
fpca_syn <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='SYNC', year_lab=TRUE, norm=FALSE)
fpca_picoe <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='PICOEC', year_lab=TRUE, norm=FALSE)
fpca_nanoe <- my_fpca(ymat_ab, info_ab, 'Abondance', sp='NANOEC', year_lab=TRUE, norm=FALSE)
fpca_other <- my_fpca(ymat_ab, info_ab, 'Abondance', sp=c('SYNC', 'PICOEC', 'NANOEC'), year_lab=TRUE)
fpca_picoe$plot
fpca_cry$plot

###export data for FDA###
export_fpca_ab <- rbind(fpca_cry$data, fpca_syn$data, fpca_picoe$data,
                        fpca_nanoe$data, fpca_pro$data)
##check inertia >5% 
ab_inertia <- data.frame('GROUPE'=c(rep('Cryptophytes', 3),
                                    rep('Synechococcus', 3),
                                    rep('Pico-eucaryotes', 3),
                                    rep('Nano-eucaryotes', 3),
                                    rep('Prochlorococcus', 3)),
                         'PC'=rep(c('PC1', 'PC2', 'PC3'), 5),
                         'INERTIA'=c(fpca_cry$inertia, fpca_syn$inertia, fpca_picoe$inertia,
                                     fpca_nanoe$inertia, fpca_pro$inertia))
ab_inertia <- ab_inertia %>% unite('VAR', c(PC, GROUPE), sep="_")
ind <- which(ab_inertia$INERTIA<5)
out <- ab_inertia$VAR[ind]

###format data for export
export_fpca_ab <- export_fpca_ab %>% 
  pivot_wider(names_from = GROUPE,
              values_from = c(PC1, PC2, PC3))

#remove inertia <5%
export_fpca_ab <- export_fpca_ab %>% dplyr::select(-out)
  write.csv(export_fpca_ab, "results/PC_AB.csv", row.names=FALSE)

###fpca on nutrients###
out1 <- which(info_nutri$GROUPE=='NH4' & 
                info_nutri$SITE==11 & 
                info_nutri$ANNEE==2018)
out2 <- which(info_nutri$GROUPE=='PO4' & 
                info_nutri$SITE==11 & 
                info_nutri$ANNEE==2014)
outs <- c(out1, out2, out3)
ymat_nutri <- ymat_nutri[,-outs]
info_nutri <- info_nutri[-outs,]
fpca_nutri <- my_fpca(ymat_nutri, info_nutri, 'Concentration', nutri=TRUE)
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

fpca_NH4$plot
fpca_NO3$plot
fpca_P$plot

###export data for FDA###
export_fpca_nutri <- rbind(fpca_NH4$data, fpca_NO3$data, fpca_NO2$data,
                           fpca_S$data, fpca_P$data)
##check inertia >5% 
nutri_inertia <- data.frame('GROUPE'=c(rep('Amonium', 3),
                                    rep('Nitrate', 3),
                                    rep('Nitrite', 3),
                                    rep('Silice', 3),
                                    rep('Phosphate', 3)),
                         'PC'=rep(c('PC1', 'PC2', 'PC3'), 5),
                         'INERTIA'=c(fpca_NH4$inertia, fpca_NO3$inertia, fpca_NO2$inertia,
                                     fpca_S$inertia, fpca_P$inertia))
nutri_inertia <- nutri_inertia %>% unite('VAR', c(PC, GROUPE), sep="_")
ind <- which(nutri_inertia$INERTIA<5)
out <- nutri_inertia$VAR[ind]

###format data for export
export_fpca_nutri <- export_fpca_nutri %>% 
  pivot_wider(names_from = GROUPE,
              values_from = c(PC1, PC2, PC3))

#remove inertia <5%
export_fpca_nutri <- export_fpca_nutri %>% dplyr::select(-out)
write.csv(export_fpca_nutri, "results/PC_NUTRI.csv", row.names=FALSE)


###fpca on diffusion###
fpca_diff <- my_fpca(ymat_diff, info_diff, 'Diffusion')
fpca_dcry <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='CRYSSC', year_lab=TRUE, norm=FALSE)
fpca_dpro <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='PROSSC', year_lab=TRUE, norm=FALSE)
fpca_dsyn <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='SYNSSC', year_lab=TRUE, norm=FALSE)
fpca_dpicoe <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='PICOESSC', year_lab=TRUE, norm=FALSE)
fpca_dnanoe <- my_fpca(ymat_diff, info_diff, 'Diffusion', sp='NANOESSC', year_lab=TRUE, norm=FALSE)
fpca_dnanoe$plot

###export data for FDA###
export_fpca_diff <- rbind(fpca_dcry$data, fpca_dsyn$data, fpca_dpicoe$data,
                        fpca_dnanoe$data, fpca_dpro$data)

##check inertia >5% 
diff_inertia <- data.frame('GROUPE'=c(rep('Cryptophytes', 3),
                                    rep('Synechococcus', 3),
                                    rep('Pico-eucaryotes', 3),
                                    rep('Nano-eucaryotes', 3),
                                    rep('Prochlorococcus', 3)),
                         'PC'=rep(c('PC1', 'PC2', 'PC3'), 5),
                         'INERTIA'=c(fpca_dcry$inertia, fpca_dsyn$inertia, fpca_dpicoe$inertia,
                                     fpca_dnanoe$inertia, fpca_dpro$inertia))
diff_inertia <- diff_inertia %>% unite('VAR', c(PC, GROUPE), sep="_")
ind <- which(diff_inertia$INERTIA<=6)
out <- diff_inertia$VAR[ind]

###format data for export
export_fpca_diff <- export_fpca_diff %>% 
  pivot_wider(names_from = GROUPE,
              values_from = c(PC1, PC2, PC3))

#remove inertia <5%
export_fpca_diff <- export_fpca_diff %>% dplyr::select(-out)

write.csv(export_fpca_diff, "results/PC_DIFF.csv", row.names=FALSE)

