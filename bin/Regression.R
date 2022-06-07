#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Regression with Confidence Intervals
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(latex2exp)
library(ggpubr)

###source code###
source("src/my_palette.R")
source("src/reg_pol.R")
source("src/weight_gauss.R")
source("src/reg_pol_taylor.R")
source("src/reg_pol_var_cond.R")

###function to compute regression and CI###

#site as ID_SITE, group as value and variable as column name
reg_CI <- function(data, site, group=NULL, variable, h, pilot_h, start='1995-01-01', ctd=FALSE){
  p <- 2 #regression order
  a <- 1 #taylor approx order
  if(is.null(group)==TRUE & ctd==FALSE){
    df <- data %>% filter(ID_SITE==site & PROF_TEXT=='S')
  }
  if(is.null(group)==FALSE){
    df <- data %>% filter(ID_SITE==site & GROUPE==group & PROF_TEXT=='S')
  }
  if(ctd){df <- data}
  df <- df %>% dplyr::select("DATE", variable)
  df <- na.omit(df)
  #date to time vector
  start <- ymd(start)
  df$DATE <- ymd(df$DATE)
  df <- df %>% filter(DATE>start)
  T0 <- df$DATE[1]
  xech <- as.numeric(df$DATE-T0)
  #values to log to avoid negative approx
  if(ctd==FALSE){
  yech <- df[,2]
  yech <- log(yech+1)
  }
  if(ctd){yech <- df[,variable]}
  #polynomial regression
  my_reg <- reg_pol(xech, yech, h, p)
  reg_time <- my_reg$X
  reg_val <- my_reg$Mp #log of value+1
  reg_date <- as.Date(reg_time, origin=T0)
  
  ##Confidence Interval
  #taylor approximation
  taylor <- reg_pol_taylor(xech, yech, pilot_h, p, a)
  ##RSME
  sigma <- taylor$sigma
  ##conditional variance matrices
  Var_cond <- var_cond(xech, yech, h, p, a, sigma)
  bhat <- taylor$biais[1,]
  Vhat <- Var_cond[1,1,]
  #student quantile
  alpha <- 0.05
  prob <- 1-alpha/2
  ld <- floor(taylor$dn) #liberty degrees
  zalpha <- qt(prob, ld) #quantile student
  #CI
  CI_up <- reg_val - bhat + zalpha * sqrt(Vhat)
  CI_low <- reg_val - bhat - zalpha * sqrt(Vhat)
  result <- data.frame("DATE" = reg_date, 
                       "TIME" = reg_time, 
                       "VALUE"= reg_val, 
                       "CI_up"= CI_up, 
                       "CI_low"= CI_low)
  colnames(result)[3] <- variable
  return(result)
}

###function to plot regression w/ CI###

#site, variable and group as display on plot
plot_reg_CI <- function(data, site, group=NULL, variable, my_points, var_points='ABONDANCE'){
  if(is.null(group)==FALSE){
  my_points$DATE <- ymd(my_points$DATE)
  my_points$ID_SITE <- factor(my_points$ID_SITE,
                              levels=c(10, 11, 12),
                              labels=c('Banyuls', 'Marseille', 'Villefranche'))
  
    my_points$GROUPE <- factor(my_points$GROUPE,
                             labels=c('Cryptophytes', 'HNA', 'LNA', 'Nano-eucaryotes',
                             'Pico-eucaryotes', 'Prochlorococcus', 'Synechococcus', 'BAC'))
    my_points <- my_points %>% filter(ID_SITE==site & GROUPE==group)
  }
  data$DATE <- ymd(data$DATE)
  colnames(data)[3] <- "VARIABLE"
  #turn log value to original value
  data$VARIABLE <- exp(data$VARIABLE)-1
  neg_reg <- which(data$VARIABLE < 0)
  data$VARIABLE[neg_reg] <- 0
  data$CI_low <- exp(data$CI_low)-1
  neg <- which(data$CI_low < 0)
  data$CI_low[neg] <- 0
  data$CI_up <- exp(data$CI_up)-1
  neg_up <- which(data$CI_up < 0)
  data$CI_up[neg_up] <- 0
  my_plot <- data %>% ggplot() + 
    geom_ribbon(aes(DATE, ymin=CI_low, ymax=CI_up), alpha=.6, fill='grey') +
    geom_line(aes(DATE, VARIABLE), col=my_palette[2], size=.8) +
    theme_light() + ylab(variable) + xlab('Temps') + 
    ggtitle(paste(site, group, sep=" - ")) +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          title=element_text(size=16))
  if(is.null(group)==FALSE){
    my_plot <- my_plot +
      geom_point(data=my_points, aes_string('DATE', var_points), size=.9) +
      geom_line(aes(DATE, VARIABLE), col=my_palette[2], size=.8)
  }
  return(my_plot)
}

###load data###
data_piconano <- read.csv("data/PICONANO_AB.csv")
data_piconano_diff <- read.csv("data/PICONANO_DIFF.csv")
data_hydro <- read.csv("data/HYDRO.csv")

###compute regression###
###ABUNDANCE
##BANYULS
#export <- reg_CI(data_piconano, site=10, group='CRYC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/cry_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='SYNC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/syn_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='PROC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/pro_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/picoe_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/nanoe_ab_banyuls.csv", row.names=FALSE)

##MARSEILLE
#export <- reg_CI(data_piconano, site=11, group='CRYC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/cry_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='SYNC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/syn_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='PROC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/pro_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/picoe_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/nanoe_ab_marseille.csv", row.names=FALSE)

##VILLEFRANCHE
#export <- reg_CI(data_piconano, site=12, group='CRYC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/cry_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='SYNC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/syn_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='PROC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/pro_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/picoe_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/nanoe_ab_villefranche.csv", row.names=FALSE)

###DIFFUSION
##BANYULS
#export <- reg_CI(data_piconano, site=10, group='CRYSSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/cry_diff_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='SYNSSC', variable='DIFFUSION', h=40, pilot_h=24)
#write.csv(export, "results/syn_diff_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='PROSSC', variable='DIFFUSION', h=35, pilot_h=30)
#write.csv(export, "results/pro_diff_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='PICOESSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/picoe_diff_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=10, group='NANOESSC', variable='DIFFUSION', h=40, pilot_h=25)
#write.csv(export, "results/nanoe_diff_banyuls.csv", row.names=FALSE)

##MARSEILLE
#export <- reg_CI(data_piconano, site=11, group='CRYSSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/cry_diff_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='SYNSSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/syn_diff_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='PROSSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/pro_diff_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='PICOESSC', variable='DIFFUSION', h=40, pilot_h=23)
#write.csv(export, "results/picoe_diff_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=11, group='NANOESSC', variable='DIFFUSION', h=40, pilot_h=18)
#write.csv(export, "results/nanoe_diff_marseille.csv", row.names=FALSE)

##VILLEFRANCHE
#export <- reg_CI(data_piconano, site=12, group='CRYSSC', variable='DIFFUSION', h=40, pilot_h=19.5)
#write.csv(export, "results/cry_diff_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='SYNSSC', variable='DIFFUSION', h=40, pilot_h=19.5)
#write.csv(export, "results/syn_diff_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='PROSSC', variable='DIFFUSION', h=30, pilot_h=22)
#write.csv(export, "results/pro_diff_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='PICOESSC', variable='DIFFUSION', h=40, pilot_h=19.5)
#write.csv(export, "results/picoe_diff_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_piconano, site=12, group='NANOESSC', variable='DIFFUSION', h=40, pilot_h=19.5)
#write.csv(export, "results/nanoe_diff_villefranche.csv", row.names=FALSE)


###NUTRIENTS
##BANYULS
#export <- reg_CI(data_hydro, site=10, variable='NH4', h=40, pilot_h=21, start='2011-09-01')
#write.csv(export, "results/nh4_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=10, variable='NO3', h=40, pilot_h=21, start='2011-11-22')
#write.csv(export, "results/no3_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=10, variable='NO2', h=40, pilot_h=18, start='2011-11-22')
#write.csv(export, "results/no2_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=10, variable='PO4', h=33, pilot_h=19, start='2011-11-22')
#write.csv(export, "results/po4_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=10, variable='SIOH4', h=33, pilot_h=19, start='2011-11-22')
#write.csv(export, "results/sioh4_banyuls.csv", row.names=FALSE)

##MARSEILLE
#export <- reg_CI(data_hydro, site=11, variable='NH4', h=45, pilot_h=21, start='2011-09-01')
#write.csv(export, "results/nh4_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=11, variable='NO3', h=45, pilot_h=20, start='2011-11-22')
#write.csv(export, "results/no3_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=11, variable='NO2', h=40, pilot_h=21, start='2011-11-22')
#write.csv(export, "results/no2_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=11, variable='PO4', h=40, pilot_h=18.5, start='2011-11-22')
#write.csv(export, "results/po4_marseille.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=11, variable='SIOH4', h=40, pilot_h=21, start='2011-11-22')
#write.csv(export, "results/sioh4_marseille.csv", row.names=FALSE)

##VILLEFRANCHE
#export <- reg_CI(data_hydro, site=12, variable='NH4', h=45, pilot_h=18, start='2011-09-01')
#write.csv(export, "results/nh4_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=12, variable='NO3', h=40, pilot_h=18, start='2011-11-22')
#write.csv(export, "results/no3_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=12, variable='NO2', h=30, pilot_h=17, start='2011-11-22')
#write.csv(export, "results/no2_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=12, variable='PO4', h=35, pilot_h=18, start='2011-11-22')
#write.csv(export, "results/po4_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data_hydro, site=12, variable='SIOH4', h=40, pilot_h=21, start='2011-11-22')
#write.csv(export, "results/sioh4_villefranche.csv", row.names=FALSE)
                                                                                                                

###BIOMASS
##BANYULS
#export <- reg_CI(data_hydro, site=10, variable='CHLA', h=40, pilot_h=21, start="2007-11-20")
#write.csv(export, "results/chla_banyuls.csv", row.names=FALSE)

##MARSEILLE
#export <- reg_CI(data_hydro, site=11, variable='CHLA', h=40, pilot_h=21, start="2007-11-20")
#write.csv(export, "results/chla_marseille.csv", row.names=FALSE)


##VILLEFRANCHE
#export <- reg_CI(data_hydro, site=12, variable='CHLA', h=35, pilot_h=21, start="2007-11-20")
#write.csv(export, "results/chla_villefranche.csv", row.names=FALSE)


###plot regression###
#abundance 
data <- read.csv("results/TS_CRY_AB_BANYULS.csv")
cry_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Cryptophytes', 'Abondance', data_piconano)
data <- read.csv("results/TS_SYN_AB_BANYULS.csv")
syn_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Synechococcus', 'Abondance', data_piconano)
data <- read.csv("results/TS_PRO_AB_BANYULS.csv")
pro_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Prochlorococcus', 'Abondance', data_piconano) + 
  scale_y_continuous(limits=c(-1, 3e4))
data <- read.csv("results/TS_PICOE_AB_BANYULS.csv")
picoe_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Pico-eucaryotes', 'Abondance', data_piconano) + 
  scale_y_continuous(limits=c(-1, 3e4))
data <- read.csv("results/TS_NANOE_AB_BANYULS.csv")
nanoe_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Nano-eucaryotes', 'Abondance', data_piconano) +  
  scale_y_continuous(limits=c(-1, 6e3))

data <- read.csv("results/TS_CRY_AB_MARSEILLE.csv")
cry_ab_mar <- plot_reg_CI(data, 'Marseille', 'Cryptophytes', 'Abondance', data_piconano) + 
  scale_y_continuous(limits=c(-1, 1.2e3))
data <- read.csv("results/TS_SYN_AB_MARSEILLE.csv")
syn_ab_mar <- plot_reg_CI(data, 'Marseille', 'Synechococcus', 'Abondance', data_piconano)
data <- read.csv("results/TS_PRO_AB_MARSEILLE.csv")
pro_ab_mar <- plot_reg_CI(data, 'Marseille', 'Prochlorococcus', 'Abondance', data_piconano)
data <- read.csv("results/TS_PICOE_AB_MARSEILLE.csv")
picoe_ab_mar <- plot_reg_CI(data, 'Marseille', 'Pico-eucaryotes', 'Abondance', data_piconano)
data <- read.csv("results/TS_NANOE_AB_MARSEILLE.csv")
nanoe_ab_mar <- plot_reg_CI(data, 'Marseille', 'Nano-eucaryotes', 'Abondance', data_piconano) 

data <- read.csv("results/TS_CRY_AB_VILLEFRANCHE.csv")
cry_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Cryptophytes', 'Abondance', data_piconano)
data <- read.csv("results/TS_SYN_AB_VILLEFRANCHE.csv")
syn_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Synechococcus', 'Abondance', data_piconano)
data <- read.csv("results/TS_PRO_AB_VILLEFRANCHE.csv")
pro_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Prochlorococcus', 'Abondance', data_piconano)
data <- read.csv("results/TS_PICOE_AB_VILLEFRANCHE.csv")
picoe_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Pico-eucaryotes', 'Abondance', data_piconano)
data <- read.csv("results/TS_NANOE_AB_VILLEFRANCHE.csv")
nanoe_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Nano-eucaryotes', 'Abondance', data_piconano) 

#diffusion
data <- read.csv("results/TS_CRY_DIFF_BANYULS.csv")
cry_diff_ban <- plot_reg_CI(data, 'Banyuls', 'Cryptophytes', 'Diffusion lumineuse') + 
  scale_y_continuous(limits=c(0, 150))
data <- read.csv("results/TS_SYN_DIFF_BANYULS.csv")
syn_diff_ban <- plot_reg_CI(data, 'Banyuls', 'Synechococcus', 'Diffusion lumineuse')
data <- read.csv("results/TS_PRO_DIFF_BANYULS.csv")
pro_diff_ban <- plot_reg_CI(data, 'Banyuls', 'Prochlorococcus', 'Diffusion lumineuse') + 
  scale_y_continuous(limits=c(0, .25))
data <- read.csv("results/TS_PICOE_DIFF_BANYULS.csv")
picoe_diff_ban <- plot_reg_CI(data, 'Banyuls', 'Pico-eucaryotes', 'Diffusion lumineuse') 
data <- read.csv("results/TS_NANOE_DIFF_BANYULS.csv")
nanoe_diff_ban <- plot_reg_CI(data, 'Banyuls', 'Nano-eucaryotes', 'Diffusion lumineuse')  
  scale_y_continuous(limits=c(-1, 6e3))

data <- read.csv("results/TS_CRY_DIFF_MARSEILLE.csv")
cry_diff_mar <- plot_reg_CI(data, 'Marseille', 'Cryptophytes', 'Diffusion lumineuse') 
data <- read.csv("results/TS_SYN_DIFF_MARSEILLE.csv")
syn_diff_mar <- plot_reg_CI(data, 'Marseille', 'Synechococcus', 'Diffusion lumineuse')
data <- read.csv("results/TS_PRO_DIFF_MARSEILLE.csv")
pro_diff_mar <- plot_reg_CI(data, 'Marseille', 'Prochlorococcus', 'Diffusion lumineuse')
data <- read.csv("results/TS_PICOE_DIFF_MARSEILLE.csv")
picoe_diff_mar <- plot_reg_CI(data, 'Marseille', 'Pico-eucaryotes', 'Diffusion lumineuse')
data <- read.csv("results/TS_NANOE_DIFF_MARSEILLE.csv")
nanoe_diff_mar <- plot_reg_CI(data, 'Marseille', 'Nano-eucaryotes', 'Diffusion lumineuse') 

data <- read.csv("results/TS_CRY_DIFF_VILLEFRANCHE.csv")
cry_diff_vil <- plot_reg_CI(data, 'Villefranche', 'Cryptophytes', 'Diffusion lumineuse')
data <- read.csv("results/TS_SYN_DIFF_VILLEFRANCHE.csv")
syn_diff_vil <- plot_reg_CI(data, 'Villefranche', 'Synechococcus', 'Diffusion lumineuse')
data <- read.csv("results/TS_PRO_DIFF_VILLEFRANCHE.csv")
pro_diff_vil <- plot_reg_CI(data, 'Villefranche', 'Prochlorococcus', 'Diffusion lumineuse')
data <- read.csv("results/TS_PICOE_DIFF_VILLEFRANCHE.csv")
picoe_diff_vil <- plot_reg_CI(data, 'Villefranche', 'Pico-eucaryotes', 'Diffusion lumineuse')
data <- read.csv("results/TS_NANOE_DIFF_VILLEFRANCHE.csv")
nanoe_diff_vil <- plot_reg_CI(data, 'Villefranche', 'Nano-eucaryotes', 'Diffusion lumineuse') 


#nutrients
data <- read.csv("results/TS_NH4_BANYULS.csv")
nh4_ban <- plot_reg_CI(data, 'Banyuls', TeX('$NH_4$'), TeX('Concentration $(\\mu M)$')) + 
  scale_y_continuous(limits=c(-0.01, 0.3))
data <- read.csv("results/TS_NO3_BANYULS.csv")
no3_ban <- plot_reg_CI(data, 'Banyuls', TeX('$NO_3$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_NO2_BANYULS.csv")
no2_ban <- plot_reg_CI(data, 'Banyuls', TeX('$NO_2$'), TeX('Concentration $(\\mu M)$')) 
data <- read.csv("results/TS_PO4_BANYULS.csv")
po4_ban <- plot_reg_CI(data, 'Banyuls', TeX('$PO_4$'), TeX('Concentration $(\\mu M)$')) + 
  scale_y_continuous(limits=c(-0.01, 0.3))
data <- read.csv("results/TS_SIOH4_BANYULS.csv")
sioh4_ban <- plot_reg_CI(data, 'Banyuls', TeX('$Si(OH)_4$'), TeX('Concentration $(\\mu M)$')) + 
  scale_y_continuous(limits=c(-0.01, 12))

data <- read.csv("results/TS_NH4_MARSEILLE.csv")
nh4_mar <- plot_reg_CI(data, 'Marseille', TeX('$NH_4$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_NO3_MARSEILLE.csv")
no3_mar <- plot_reg_CI(data, 'Marseille', TeX('$NO_3$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_NO2_MARSEILLE.csv")
no2_mar <- plot_reg_CI(data, 'Marseille', TeX('$NO_2$'), TeX('Concentration $(\\mu M)$')) 
data <- read.csv("results/TS_PO4_MARSEILLE.csv")
po4_mar <- plot_reg_CI(data, 'Marseille', TeX('$PO_4$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_SIOH4_MARSEILLE.csv")
sioh4_mar <- plot_reg_CI(data, 'Marseille', TeX('$Si(OH)_4$'), TeX('Concentration $(\\mu M)$')) 

data <- read.csv("results/TS_NH4_VILLEFRANCHE.csv")
nh4_vil <- plot_reg_CI(data, 'Villefranche', TeX('$NH_4$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_NO3_VILLEFRANCHE.csv")
no3_vil <- plot_reg_CI(data, 'Villefranche', TeX('$NO_3$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_NO2_VILLEFRANCHE.csv")
no2_vil <- plot_reg_CI(data, 'Villefranche', TeX('$NO_2$'), TeX('Concentration $(\\mu M)$')) 
data <- read.csv("results/TS_PO4_VILLEFRANCHE.csv")
po4_vil <- plot_reg_CI(data, 'Villefranche', TeX('$PO_4$'), TeX('Concentration $(\\mu M)$'))
data <- read.csv("results/TS_SIOH4_VILLEFRANCHE.csv")
sioh4_vil <- plot_reg_CI(data, 'Villefranche', TeX('$Si(OH)_4$'), TeX('Concentration $(\\mu M)$')) 

##biomass
data <- read.csv("results/TS_CHLA_BANYULS.csv")
chla_ban <- plot_reg_CI(data, 'Banyuls', 'Chlorophylle a', TeX('Concentration $(\\mu g/L)$'))
data <- read.csv("results/TS_CHLA_MARSEILLE.csv")
chla_mar <- plot_reg_CI(data, 'Marseille', 'Chlorophylle a', TeX('Concentration $(\\mu g/L)$'))
data <- read.csv("results/chla_villefranche.csv")
chla_vil <- plot_reg_CI(data, 'Villefranche', 'Chlorophylle a', TeX('Concentration $(\\mu g/L)$'))

ggarrange(syn_ab_ban + xlab('') + scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)), 
          syn_ab_mar + ylab('') + xlab('') + scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)),
          syn_ab_vil + ylab('') + xlab('')+ scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)),
          cry_ab_ban + xlab('') + scale_y_continuous(limits=c(0, 1250)),
          cry_ab_mar + ylab('') + xlab('')+ scale_y_continuous(limits=c(0, 1250)),
          cry_ab_vil  + xlab('')+ ylab('')+ scale_y_continuous(limits=c(0, 1250)),
          pro_ab_ban + xlab('')+ scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          pro_ab_mar + ylab('') + xlab('')+ scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          pro_ab_vil +  xlab('')+ ylab('')+ scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          nanoe_ab_ban + xlab('')+ scale_y_continuous(limits=c(0, 7000)),
          nanoe_ab_mar + xlab('') + ylab('')+ scale_y_continuous(limits=c(0, 7000)),
          nanoe_ab_vil  + xlab('') + ylab('')+ scale_y_continuous(limits=c(0, 7000)),
          picoe_ab_ban + scale_y_continuous(limits=c(0, 35000)),
          picoe_ab_mar + ylab('')+ scale_y_continuous(limits=c(0, 35000)),
          picoe_ab_vil + ylab('')+ scale_y_continuous(limits=c(0, 35000)),
          ncol=3, nrow=5)

ggarrange(syn_ab_ban + xlab('') + scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)), 
          syn_ab_mar + ylab('') + xlab('') + scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)),
          syn_ab_vil + ylab('') + xlab('')+ scale_y_continuous(limits=c(0, 120000), breaks=seq(0, 120000, length.out=4)),
          cry_ab_ban + xlab('') + scale_y_continuous(limits=c(0, 1250)),
          cry_ab_mar + ylab('') + xlab('')+ scale_y_continuous(limits=c(0, 1250)),
          cry_ab_vil  + xlab('')+ ylab('')+ scale_y_continuous(limits=c(0, 1250)),
          pro_ab_ban +  scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          pro_ab_mar + ylab('') +  scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          pro_ab_vil +   ylab('')+ scale_y_continuous(limits=c(0, 45000), breaks=seq(0, 45000, length.out=4)),
          ncol=3, nrow=3)
          

###PC from CTD FPCA###
data_CTD <- read.csv('results/DATA_PC_CTD_MARSEILLE.csv')
#ctd_mar <- reg_CI(data_CTD, site=11, variable='PC1', h=40, pilot_h=19.5, ctd=TRUE)
#write.csv(ctd_mar, 'results/TS_PC1_CTD_MARSEILLE.csv', row.names=FALSE)
data_CTD$DATE <- ymd(data_CTD$DATE)
date_range <- ymd(c('1995-01-01', '2022-01-01'))
ctd_mar <- read.csv('results/TS_PC1_CTD_MARSEILLE.csv')
ctd_mar$DATE <- ymd(ctd_mar$DATE)
#ctd_mar2 <- reg_CI(data_CTD, site=11, variable='PC2', h=40, pilot_h=19.5, ctd=TRUE)
#write.csv(ctd_mar2, 'results/TS_PC2_CTD_MARSEILLE.csv', row.names=FALSE)
ctd_mar2 <- read.csv('results/TS_PC2_CTD_MARSEILLE.csv')
ctd_mar2$DATE <- ymd(ctd_mar2$DATE)
A <- ctd_mar %>% ggplot() + 
  geom_ribbon(aes(DATE, ymin=CI_low, ymax=CI_up), alpha=.6, fill='grey') +
  geom_point(data=data_CTD, aes(DATE, PC1), size=.9) +
  geom_line(aes(DATE, PC1), col=my_palette[2], size=.8) +
  scale_y_continuous(limits=c(-50, 50)) +
  scale_x_date(limits=date_range, 
               breaks=seq.Date(ymd("1995-01-01"), ymd("2022-01-01"), by='2 year'),
               labels=seq(1995, 2022, 2)) +
  theme_light() + ylab('PC1') + xlab('Temps') + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18))
B <- ctd_mar2 %>% ggplot() + 
  geom_ribbon(aes(DATE, ymin=CI_low, ymax=CI_up), alpha=.6, fill='grey') +
  geom_point(data=data_CTD, aes(DATE, PC1), size=.9) +
  geom_line(aes(DATE, PC1), col=my_palette[2], size=.8) +
  scale_y_continuous(limits=c(-30, 50)) +
  scale_x_date(limits=date_range, 
               breaks=seq.Date(ymd("1995-01-01"), ymd("2022-01-01"), by='2 year'),
               labels=seq(1995, 2022, 2)) +
  theme_light() + ylab('PC1') + xlab('Temps') + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18))
