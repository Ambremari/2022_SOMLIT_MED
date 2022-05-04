#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#FDA
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(MASS)
library(missMDA)
library(ggforce)

###source code###
source("src/my_palette.R")

###load data###
data_ab <- read.csv('results/PC_AB.csv')
data_nutri <- read.csv('results/PC_NUTRI.csv')
data_diff <- read.csv('results/PC_DIFF.csv')

###fda plots function###
my_fda <- function(data, variable, nutri=FALSE){
  mat <- data[,-(1:2)]
  ##NA treatment
  imputed <- imputePCA(mat, ncp=2)
  #fda
  input <- data.frame("SITE"=data[,2], imputed$completeObs)
  var <- lda(SITE~., input)
  ind <- predict(var, data=input)
  p_axe1 <- round(prop.table(var$svd^2)[1]*100, 2)
  p_axe2 <- round(prop.table(var$svd^2)[2]*100, 2)
  my_df <- data.frame("ANNEE"=data[,1],"SITE"=data[,2], ind$x)
  my_df %>% ggplot() + geom_point(aes(LD1, LD2, col=factor(SITE))) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    theme_light()
  Xacp1 <- my_df %>% filter(SITE=="Banyuls")
  Xacp2 <- my_df %>% filter(SITE=="Marseille")
  Xacp3 <- my_df %>% filter(SITE=="Villefranche")
  g1 <- apply(Xacp1[,3:4], 2, mean)
  g2 <- apply(Xacp2[,3:4], 2, mean)
  g3 <- apply(Xacp3[,3:4], 2, mean)
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
              fontface='bold', size=4.5, vjust=1, hjust=1) +
    geom_segment(data=ban, aes(g1[1], g1[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    geom_segment(data=mar, aes(g2[1], g2[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    geom_segment(data=vil, aes(g3[1], g3[2], xend=LD1, yend=LD2, col=factor(SITE))) +
    scale_color_manual(values=my_palette[2:4], name='Site') +
    xlab(paste('Axe 1 ', p_axe1, '%')) + ylab(paste('Axe 2 ', p_axe2, '%')) + ggtitle(variable)+
    theme(legend.text = element_text(size=13),
          legend.position='none',
          legend.title = element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  coord <- data.frame('VAR'=row.names(var$scaling), var$scaling)
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
    scale_fill_manual(values=my_palette, name=leg_name)
  return(list(acp=acp_plot, var=var_plot, data=my_df))
}

####ABUNDANCE
my_fda(data_ab, 'Abondance')$acp
x11()
my_fda(data_ab, 'Abondance')$var

####NUTRIENTS
x11()
my_fda(data_nutri, 'Nutriments', nutri=TRUE)$acp
x11()
my_fda(data_nutri, 'Nutriments', nutri=TRUE)$var

####DIFFUSION
x11()
my_fda(data_diff, 'Diffusion', nutri=TRUE)$acp
x11()
my_fda(data_diff, 'Diffusion', nutri=TRUE)$var

