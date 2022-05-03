#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#FDA
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(FactoMineR)
library(missMDA)
library(ggforce)

###source code###
source("src/my_palette.R")

###load data###
data_ab <- read.csv('results/PC_AB.csv')
data_nutri <- read.csv('results/PC_NUTRI.csv')

###fda plots function###
my_fda <- function(data, variable, nutri=FALSE){
  mat <- data[,-(1:2)]
  ##NA treatment
  imputed <- imputePCA(mat, ncp=2)
  acp <- PCA(imputed$completeObs, ncp=2, graph=FALSE)
  p_axe1 <- round(acp$eig[1,2], 2)
  p_axe2 <- round(acp$eig[2,2], 2)
  my_df <- data.frame("SITE"=data[,2], acp$ind$coord)
  my_df %>% ggplot() + geom_point(aes(Dim.1, Dim.2, col=factor(SITE))) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    theme_light()
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
  acp_plot <- my_df %>% ggplot() + geom_point(aes(Dim.1, Dim.2, col=factor(SITE))) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    theme_light() +
    geom_point(data=grav_center, aes(Dim.1, Dim.2, col=factor(SITE)), shape=3) +
    geom_segment(data=ban, aes(g1[1], g1[2], xend=Dim.1, yend=Dim.2, col=factor(SITE))) +
    geom_segment(data=mar, aes(g2[1], g2[2], xend=Dim.1, yend=Dim.2, col=factor(SITE))) +
    geom_segment(data=vil, aes(g3[1], g3[2], xend=Dim.1, yend=Dim.2, col=factor(SITE))) +
    scale_color_manual(values=my_palette[2:4], name='Site') +
    xlab(paste('Axe 1 ', p_axe1, '%')) + ylab(paste('Axe 2 ', p_axe2, '%')) + ggtitle(variable)+
    theme(legend.text = element_text(size=13),
          legend.title = element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  coord <- data.frame('VAR'=row.names(acp$var$coord), acp$var$coord)
  coord <- coord %>% separate('VAR', c('PC', 'GROUPE'), '_')
  contrib <-  data.frame('VAR'=row.names(acp$var$coord), acp$var$contrib) 
  contrib <- contrib %>% separate('VAR', c('PC', 'GROUPE'), '_')
  leg_name <- 'Groupe'
  if(nutri==TRUE){
    contrib$GROUPE <- factor(contrib$GROUPE, 
                             labels=c('NH4', 'NO3', 'NO2', 'PO4', 'SiOH(4)'))
    leg_name <- 'Nutriment'
  }
  var_plot <- coord %>% ggplot() + 
    geom_segment(aes(0, 0, xend=Dim.1, yend=Dim.2)) +
    geom_vline(xintercept=0, linetype='dashed') +
    geom_hline(yintercept=0, linetype='dashed') +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), inherit.aes = FALSE) +
    geom_label(aes(Dim.1, Dim.2, label=PC, fill=GROUPE), 
               color='white', size=4.5) +
    theme_light() + xlab(paste('Axe 1 ', p_axe1, '%')) + ylab(paste('Axe 2 ', p_axe2, '%')) + ggtitle(variable) +
    scale_fill_manual(values=my_palette, name=leg_name)
  contrib <- contrib %>% unite('VAR', c('PC', 'GROUPE'), sep='_', remove=FALSE)
  axe_un <- contrib %>% arrange(desc(Dim.1)) %>% 
    mutate(VAR=factor(VAR, levels=VAR)) %>%
    ggplot() + theme_minimal() + xlab('Variable') +
    geom_bar(aes(VAR, Dim.1, fill=GROUPE), 
             stat='identity', width=.6, alpha=.8) +
    scale_fill_manual(values=my_palette, name=leg_name) +
    theme(axis.text.x = element_text(size=9, angle=45, hjust=1),
          legend.text = element_text(size=13),
          legend.title = element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  axe_deux <- contrib %>% arrange(desc(Dim.2)) %>% 
    mutate(VAR=factor(VAR, levels=VAR)) %>%
    ggplot() + theme_minimal() + xlab('Variable') +
    geom_bar(aes(VAR, Dim.2, fill=GROUPE), 
             stat='identity', width=.6, alpha=.8) +
    scale_fill_manual(values=my_palette, name=leg_name) +
    theme(axis.text.x = element_text(size=9, angle=45, hjust=1),
          legend.text = element_text(size=13),
          legend.title = element_text(size=14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
  if(nutri==FALSE){
    axe_un <- axe_un + scale_x_discrete(labels=abbreviate)
    axe_deux <- axe_deux + scale_x_discrete(labels=abbreviate) 
  }
  contrib_plot <- ggarrange(NULL, axe_un, NULL, axe_deux, 
                            ncol=1, nrow=4,
                            heights=c(0.1, 1, 0.1, 1),
                            labels=c('', 'Axe 1', '', 'Axe 2'),
                            font.label=list(size=12, face='plain'),
                            common.legend=TRUE, legend='bottom',
                            vjust=0.6)
  return(list(acp=acp_plot, var=var_plot, contrib=contrib_plot))
}

####ABUNDANCE
my_fda(data_ab, 'Abondance')$acp
x11()
my_fda(data_ab, 'Abondance')$var
x11()
my_fda(data_ab, 'Abondance')$contrib

####ANUTRIENTS
my_fda(data_nutri, 'Nutriments', nutri=TRUE)$acp
x11()
my_fda(data_nutri, 'Nutriments', nutri=TRUE)$var
x11()
my_fda(data_nutri, 'Nutriments', nutri=TRUE)$contrib