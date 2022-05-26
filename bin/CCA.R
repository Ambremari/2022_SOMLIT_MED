#################
#SOMLIT MED - 2022/05/03
#Mathilde Couteyen Carpaye 
#CCA
################


###working directory###
setwd("C:/Users/ambre/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")


###load packages###
library(tidyverse)
library(CCA)
library(ggpubr)

###source codes###
source('src/my_palette.R')


###load data###
ab <- read.csv('results/PC_AB.csv')
nutri <- read.csv('results/PC_NUTRI.csv')
diff <- read.csv('results/PC_DIFF.csv')


###CCA FUNCTION ###
my_cca <- function(data1, data2, variable, dim1, dim2){
  data <- left_join(data1, data2, by=c('ANNEE', 'SITE'))
  data1 <- data[, 3:7]
  data2 <- data[,18: 22]
  cca <- cc(data2, data1)
  cor1 <- round(cca$cor[dim1], 3)
  cor2 <- round(cca$cor[dim2], 3)
  CC_X <- as.matrix(data2) %*% cca$xcoef
  my_df <- data.frame('SITE'=data$SITE, 'CC1_x'=CC_X[,dim1], 'CC2_x'=CC_X[,dim2])
  x_var <- rbind(t(t(cca$scores$corr.X.xscores[,dim1])),
                t(t(cca$scores$corr.Y.xscores[,dim1])))
  y_var <- rbind(t(t(cca$scores$corr.X.xscores[,dim2])),
               t(t(cca$scores$corr.Y.xscores[,dim2])))
  my_var <- c('NH4', 'NO3', 'NO2',  'Si(OH)4', 'PO4', 'CRY', 'SYN', 'PICOE', 'NANOE', 'PRO')
  var_df <- data.frame('VAR'=my_var, 'X'=x_var, 'Y'=y_var)
  m <- length(cca$cor)
  cor_df <- data.frame('CC'=1:m, 'cor'=cca$cor)
  cor_plot <- ggplot() + 
    geom_bar(data=cor_df, aes(CC, cor), stat='identity') +
    theme_light() + xlab('Dimension') + ylab('Corrélation canonique')
  ind_plot <- my_df %>% ggplot() + theme_light() +
   geom_point(aes(CC1_x, CC2_x, col=factor(SITE)), size=3.5) +
   scale_color_manual(values=my_palette[2:4], name='Site') +
   geom_vline(xintercept=0, linetype='dashed') +
   geom_hline(yintercept=0, linetype='dashed') +
    ggtitle(variable) +
   xlab(paste('CC', dim1, " ", cor1)) + ylab(paste('CC', dim2, " ", cor2)) + 
    theme(legend.text = element_text(size=18),
          legend.title = element_text(size=18),
          axis.text=element_text(size=16),
          axis.title=element_text(size=18),
          legend.position='bottom',
          title=element_text(size=16))
  var_plot <- var_df %>% ggplot() + theme_light() +
   geom_point(aes(X, Y), shape=c(rep(16, 5), rep(2, 5)), size=2.5) +
   geom_text(aes(X, Y, label=VAR), col=c(rep(my_palette[2], 5), rep(my_palette[3], 5)), 
            vjust=1.5, size=6) +
   geom_segment(data=var_df[1:5,], aes(x=0, y=0, xend=X, yend=Y)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    xlab(paste('CC', dim1, " ", cor1)) + ylab(paste('CC', dim2, " ", cor2)) + 
    theme(legend.text = element_text(size=16),
          legend.title = element_text(size=16),
          axis.text=element_text(size=16),
          axis.title=element_text(size=18))
  return(list('var_plot'=var_plot, 'ind_plot'=ind_plot, 'cor_plot'=cor_plot))
}

###CCA###
cca_ab <- my_cca(ab, nutri, 'Abondance', 1, 2)
ggarrange(cca_ab$ind_plot, NULL,  cca_ab$var_plot, ncol=3, nrow=1, 
          widths=c(1, .1, 1))
cca_diff <- my_cca(diff, nutri, 'Diffusion', 1, 5)
ggarrange(cca_diff$ind_plot, NULL,  cca_diff$var_plot, ncol=3, nrow=1, 
          widths=c(1, .1, 1))
