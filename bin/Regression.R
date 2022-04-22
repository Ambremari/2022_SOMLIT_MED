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

###source code###
source("src/my_palette.R")
source("src/reg_pol.R")
source("src/weight_gauss.R")
source("src/reg_pol_taylor.R")
source("src/reg_pol_var_cond.R")

###function to compute regression and CI###

#site as ID_SITE, group as value and variable as column name
reg_CI <- function(data, site, group, variable, h, pilot_h){
  p <- 2 #regression order
  a <- 1 #taylor approx order
  if(is.null(group)==TRUE){
    df <- data %>% filter(ID_SITE==site & PROF_TEXT=='S')
  }
  if(is.null(group)==FALSE){
    df <- data %>% filter(ID_SITE==site & GROUPE==group & PROF_TEXT=='S')
    }
  df <- df %>% dplyr::select("DATE", variable)
  #date to time vector
  df$DATE <- ymd(df$DATE)
  T0 <- df$DATE[1]
  xech <- as.numeric(df$DATE-T0)
  #values to log to avoid negative approx
  y_var <- df[,2]
  yech <- log(y_var+1)
  
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
plot_reg_CI <- function(data, site, group, variable){
  data$DATE <- ymd(data$DATE)
  colnames(data)[3] <- "VARIABLE"
  #turn log value to original value
  data$VARIABLE <- exp(data$VARIABLE)-1
  my_plot <- data %>% ggplot() + 
    geom_ribbon(aes(DATE, ymin=exp(CI_low), ymax=exp(CI_up)), alpha=.6, fill='grey') +
    geom_line(aes(DATE, VARIABLE), col=my_palette[2], size=.8) +
    theme_light() + ylab(variable) + xlab('Temps') + 
    ggtitle(paste(site, group, sep=" - "))
  return(my_plot)
}

###compute regression###
data <- read.csv("data/PICONANO_AB.csv")

##BANYULS
#export <- reg_CI(data, site=10, group='CRYC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/cry_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data, site=10, group='SYNC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/syn_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data, site=10, group='PROC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/pro_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data, site=10, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/picoe_ab_banyuls.csv", row.names=FALSE)
#export <- reg_CI(data, site=10, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=25)
#write.csv(export, "results/nanoe_ab_banyuls.csv", row.names=FALSE)

##MARSEILLE
#export <- reg_CI(data, site=11, group='CRYC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/cry_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data, site=11, group='SYNC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/syn_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data, site=11, group='PROC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/pro_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data, site=11, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/picoe_ab_marseille.csv", row.names=FALSE)
#export <- reg_CI(data, site=11, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=23)
#write.csv(export, "results/nanoe_ab_marseille.csv", row.names=FALSE)

##VILLEFRANCHE
#export <- reg_CI(data, site=12, group='CRYC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/cry_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data, site=12, group='SYNC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/syn_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data, site=12, group='PROC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/pro_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data, site=12, group='PICOEC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/picoe_ab_villefranche.csv", row.names=FALSE)
#export <- reg_CI(data, site=12, group='NANOEC', variable='ABONDANCE', h=40, pilot_h=19.5)
#write.csv(export, "results/nanoe_ab_villefranche.csv", row.names=FALSE)


###plot regression###
data <- read.csv("results/TS_CRY_AB_BANYULS.csv")
cry_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Cryptophytes', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_SYN_AB_BANYULS.csv")
syn_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Synechococcus', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_PRO_AB_BANYULS.csv")
pro_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Prochlorococcus', 'Abondance (cellules/mL)') + 
  scale_y_continuous(limits=c(-1, 3e4))
data <- read.csv("results/TS_PICOE_AB_BANYULS.csv")
picoe_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Pico-eucaryotes', 'Abondance (cellules/mL)') + 
  scale_y_continuous(limits=c(-1, 3e4))
data <- read.csv("results/TS_NANOE_AB_BANYULS.csv")
nanoe_ab_ban <- plot_reg_CI(data, 'Banyuls', 'Nano-eucaryotes', 'Abondance (cellules/mL)') +  
  scale_y_continuous(limits=c(-1, 6e3))

data <- read.csv("results/TS_CRY_AB_MARSEILLE.csv")
cry_ab_mar <- plot_reg_CI(data, 'Marseille', 'Cryptophytes', 'Abondance (cellules/mL)') + 
  scale_y_continuous(limits=c(-1, 1.2e3))
data <- read.csv("results/TS_SYN_AB_MARSEILLE.csv")
syn_ab_mar <- plot_reg_CI(data, 'Marseille', 'Synechococcus', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_PRO_AB_MARSEILLE.csv")
pro_ab_mar <- plot_reg_CI(data, 'Marseille', 'Prochlorococcus', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_PICOE_AB_MARSEILLE.csv")
picoe_ab_mar <- plot_reg_CI(data, 'Marseille', 'Pico-eucaryotes', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_NANOE_AB_MARSEILLE.csv")
nanoe_ab_mar <- plot_reg_CI(data, 'Marseille', 'Nano-eucaryotes', 'Abondance (cellules/mL)') 

data <- read.csv("results/TS_CRY_AB_VILLEFRANCHE.csv")
cry_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Cryptophytes', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_SYN_AB_VILLEFRANCHE.csv")
syn_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Synechococcus', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_PRO_AB_VILLEFRANCHE.csv")
pro_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Prochlorococcus', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_PICOE_AB_VILLEFRANCHE.csv")
picoe_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Pico-eucaryotes', 'Abondance (cellules/mL)')
data <- read.csv("results/TS_NANOE_AB_VILLEFRANCHE.csv")
nanoe_ab_vil <- plot_reg_CI(data, 'Villefranche', 'Nano-eucaryotes', 'Abondance (cellules/mL)') 
