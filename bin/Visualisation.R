#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#First visualisation
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

###load data###
data_hydro <- read.csv("data/HYDRO.csv")
data_piconano <- read.csv("data/PICONANO.csv")
data_ab <- read.csv("data/PICONANO_AB.csv")

###format data###
data_hydro$DATE <- ymd(data_hydro$DATE)
data_piconano$DATE <- ymd(data_piconano$DATE)
data_hydro$ID_SITE <- factor(data_hydro$ID_SITE, 
                             levels=c(10, 11, 12),
                             labels=c("Banyuls", "Marseille", "Villefranche"))

###plot function###
plot_var <- function(data, col_name, variable, protocole=NULL, out=NULL, log_trans=TRUE){
  if(is.null(out)){
    out <- data.frame("DATE"=ymd('2012-01-01'), as.numeric(NA))
    colnames(out) <- c('DATE', col_name)
  }
  if(log_trans){
    id_out <- which(colnames(out)==col_name)
    id <- which(colnames(data)==col_name)
    data[,id] <- log(data[,id]+1)
    out[,id_out] <- log(out[,id_out]+1)
  }
  protocole <- ymd(protocole)
  my_plot <- data %>% ggplot(aes_string("DATE", col_name)) + 
    geom_point(alpha=.6) +
    geom_point(data=out, aes_string("DATE", col_name), col='red', alpha=.9) +
    geom_vline(xintercept = protocole, linetype='dashed', col=my_palette[3], size=.6) +
    geom_vline(xintercept = ymd('2012-01-01'), linetype='dashed', col=my_palette[2], size=.5, alpha=.7) +
    theme_light() + xlab('Temps') + ylab(variable) + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.text.x = element_text(size = 14),
          plot.title = element_text(size=16))
  return(my_plot)
}
plot_ser <- function(data, col_name, variable){
  first_day <- seq.Date(ymd("2012-01-01"), ymd("2022-01-01"), by='year')
  id <- which(colnames(data)==col_name)
  data[,id] <- log(data[,id]+1)
  my_plot <- data %>% ggplot(aes_string("DATE", col_name)) + 
    geom_segment(aes_string(x="DATE", y=0, xend="DATE", yend=col_name), col='gray48', size=1.4) +
    theme_light() + xlab('Temps') + ylab(variable) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.text.x = element_text(size = 14),
          plot.title = element_text(size=16)) +
    geom_vline(xintercept=first_day, linetype='dashed', col=my_palette[2], lwd=0.7)
}

##surface values###
hydro_surf <- data_hydro %>% filter(PROF_TEXT=='S')
profs <- unique(hydro_surf$PROF_NUM)
checkpt <- sum(profs %in% c(1,3))
if(checkpt != 2){stop("Error in surface value")}

out_sal <- hydro_surf %>% filter(S<32)
sal <- plot_var(hydro_surf, 'S', 'Salinité', out=out_sal, log_trans=FALSE)
oxy <- plot_var(hydro_surf, 'O', 'log(Oxygène dissous)', "2008-04-01")
ph <- plot_var(hydro_surf, 'PH', 'pH', "2016-10-19", log_trans=FALSE)
nh4 <- plot_var(hydro_surf, 'NH4', TeX('$log(NH_4)$'), "2011-09-01")
no3 <- plot_var(hydro_surf, 'NO3', TeX('$log(NO_3)$'), "2011-11-22")
out_no2 <- hydro_surf %>% filter(NO2>7.5)
no2 <- plot_var(hydro_surf, 'NO2', TeX('$log(NO_2)$'), "2011-11-22", out_no2)
po4 <- plot_var(hydro_surf, 'PO4', TeX('$log(PO_4)$'), "2011-11-22")
sioh4 <- plot_var(hydro_surf, 'SIOH4', TeX('$log(Si(OH)_4)$'), "2011-11-22")
cop <- plot_var(hydro_surf, 'COP', TeX('log(COP)'), "2008-07-29")
nop <- plot_var(hydro_surf, 'NOP', TeX('log(NOP)'), "2008-07-29")
mes <- plot_var(hydro_surf, 'MES', TeX('log(MES)'), "2008-02-01")
chla <- plot_var(hydro_surf, 'CHLA', TeX('log(Chlorophylle a)'), "2007-11-20")

ggarrange(sal + xlab(''), oxy+ xlab(''), ph+ xlab(''), 
          nh4+ xlab(''), no3+ xlab(''), no2+ xlab(''),
          po4+ xlab(''), sioh4+ xlab(''), cop, nop, mes, chla, 
          ncol=4, nrow=3, common.legend = TRUE)

###bottom value###
hydro_bot <- data_hydro %>% filter(PROF_TEXT=='F')
profs <- unique(hydro_bot$PROF_NUM)
ind <- which(profs !=24)
checkpt <- profs[ind] >=49
if(sum(checkpt)!=length(checkpt)){stop("Error in bottom value")}


out_tf <- hydro_bot %>% filter(ID_SITE=="Villefranche" & DATE<ymd("2003-01-01"))
out_ts <- hydro_surf %>% filter(ID_SITE=="Villefranche" & DATE<ymd("2003-01-01"))
t_f <- plot_var(hydro_bot, 'T', 'TempÃ©rature dans le fond', out=out_tf) + 
  facet_wrap(.~ID_SITE) 
t_s <- plot_var(hydro_surf, 'T', 'TempÃ©rature en surface', out=out_ts) + 
  facet_wrap(.~ID_SITE)
ggarrange(t_f, t_s, nrow=2)


###piconano###
piconano_m <- data_piconano %>% filter(ID_SITE==11)
abc <- plot_ser(piconano_m, "CRYC", "Abondance") + ggtitle("Marseille - Cryptophytes")
diffc <- plot_ser(piconano_m, "CRYSSC", "Diffusion lumineuse") + ggtitle("") 
flrc <- plot_ser(piconano_m, "CRYFLR", "Auto-fluorescence rouge") + ggtitle("") 
floc <- plot_ser(piconano_m, "CRYFLO", "Auto-fluorescence orange") + ggtitle("") 
abs <- plot_ser(piconano_m, "SYNC", "Abondance") + ggtitle("Marseille - Synechococcus") 
diffs <- plot_ser(piconano_m, "SYNSSC", "Diffusion lumineuse") + ggtitle("")
flrs <- plot_ser(piconano_m, "SYNFLR", "Auto-fluorescence rouge") + ggtitle("")
flos <- plot_ser(piconano_m, "SYNFLO", "Auto-fluorescence orange") + ggtitle("")
ggarrange(abc, diffc, flrc, floc, abs, diffs, flrs, flos, nrow=2, ncol=4)

##diffusion and fluo
phyto_diff <- c("CRYSSC", 'SYNSSC', 'PROSSC',  
                'PICOESSC', 'NANOESSC')
phyto_fluo <- c("CRYFLR", 'SYNFLR', 'PROFLR', 'PICOEFLR', 'NANOEFLR')
dat_fluo <- data_piconano %>% dplyr::select("X", phyto_fluo) %>% 
  pivot_longer(phyto_fluo, names_to='GROUPE', values_to='FLUO')
dat_diff <- data_piconano %>% dplyr::select("X", phyto_diff) %>%
  pivot_longer(phyto_diff, names_to='GROUPE', values_to='DIFF')
data_phyto <- cbind(dat_diff, "FLUO"=dat_fluo$FLUO)
data_phyto %>% ggplot() + 
  geom_point(aes(FLUO, DIFF, col=GROUPE), size=2, alpha=.6) +
  theme_light() + scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values=my_palette, name='Groupe',
                     labels=c("Cryptophytes", "Nano-eucaryotes", "Pico-eucaryotes", 
                              "Prochlorococcus", "Synechococcus")) + xlab('Auto-fluroresence rouge') +
  ylab('Diffusion lumineuse') +
  theme(legend.text = element_text(size=13),
        legend.title = element_text(size=14),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))



