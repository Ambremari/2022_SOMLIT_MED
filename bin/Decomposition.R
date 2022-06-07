#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Multi-Seasonal Trend decomposition
################

###working directory###
setwd("C:/Users/ambre/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(forecast)
library(latex2exp)
library(astsa)

###source code###
source("src/my_palette.R")

###plot periodogram function###
plot_period <- function(data, variable, site, var_name, sp=NULL, k=c(1,1), tp=0.5, K){
  data$DATE <- ymd(data$DATE)
  ystart <- year(data$DATE[1])
  dstart <- yday(data$DATE[1])
  var <- data %>% dplyr::select(variable)
  my_ts <- ts(var, start=c(ystart, dstart), frequency=365)
  my_tsc <- my_ts-mean(my_ts) #centering
  n <- length(my_tsc)
  Per <- Mod(fft(my_tsc))^2/n
  freq <- (1:n -1)/n*365
  kern <- kernel("modified.daniell", k)
  my_per <- mvspec(my_ts, kernel=kern, taper=tp, plot=FALSE)
  alpha <- 0.5/K
  bdw <- round(my_per$bandwidth, 1)
  dgf <- round(my_per$df, 3)
  Lh <- round(my_per$Lh, 1)
  U <- qchisq(alpha/2, 2*Lh)
  L <- qchisq(1-alpha/2, 2*Lh)
  fun_per <- splinefun(my_per$freq, my_per$spec)
  x <- seq(0.1,5,0.001)
  dfun <-  fun_per(x, deriv=1)
  bd <- 0.3
  Roots <- NULL
  low <- 0.11
  upp <- bd
  while(length(Roots)<K){
    root <- uniroot(splinefun(x, dfun), c(low,upp), extendInt="yes")$root
    root <- round(root, 2)
    if(root %in% Roots==FALSE & root>0.12 & root<3.5 & fun_per(root, deriv=2)<0)
      Roots <- c(Roots, root)
    low <- upp
    upp <- upp+bd
  }
  CI_up <- 2*Lh*fun_per(Roots)/U
  CI_low <- 2*Lh*fun_per(Roots)/L
  par(mar=c(5,5,5,0.5), mgp=c(3,1,0), cex.lab=2.5, cex.axis=2, cex.main=3)
  plot.new()
  grid(lty=1)
  par(new=TRUE)
  plot(freq, Per, type='o', xlab='Fréquence', 
       ylab='Périodogramme', lwd=1, xlim=c(0,5),
       main=paste(site, "-", var_name, sp))
  plot.new()
  grid(lty=1)
  par(new=TRUE)
  plot(my_per$freq, my_per$spec, type='l', xlim=c(0,5), 
       xlab='Fréquence', ylab='Périodogramme moyen', 
       main=paste("fenêtre = ", bdw))
  abline(h=CI_low, lty=2, col=2:(K+1))
  legend("topright", lty=2, legend = Roots, col=2:(K+1), cex=2)
  }


###MSTL function###
my_mstl <- function(data, variable, frequency, seasons, sw, tw){
  my_start <- decimal_date(ymd(data$DATE[1]))
  var <- data %>% dplyr::select(variable)
  my_ts <- msts(var, seasonal.periods=seasons,
                ts.frequency=frequency, start=my_start)
  fit <- mstl(my_ts, 
              s.window=sw, t.window=tw, iterate=2) 
  
  fit_df <- cbind(as.data.frame(fit), "DATE"=ymd(data$DATE))
  all_seas <- NULL
  n <- length(seasons)
  for(i in 1:n){
    seas <- paste('Seasonal', seasons[i], sep="")
    all_seas <- c(all_seas, seas)
  }
  pred <- fit_df %>% dplyr::select('Trend', all_seas)
  fit_df <- fit_df %>% mutate("Predict"=apply(pred, 1, sum))
  return(fit_df)
}

###Decorrelation plot###
plot_cor <- function(data){
  tau <- 0:25
  RHO <- NULL
  for(i in tau){
    n <- nrow(data)
    end <- n-i
    start <- 1+i
    xt <- data$Remainder[1:end] #x(t)
    xtt <- data$Remainder[start:n]
    rho <- cor(xt, xtt)
    RHO <- c(RHO, rho)
  }
  corr <- data.frame("LAG"=tau, "RHO"=RHO)
  my_plot <- corr %>% ggplot() + geom_point(aes(LAG, RHO)) +
    geom_hline(yintercept=-0.087, col='red', linetype='dashed') +
    geom_hline(yintercept=0.087, col='red', linetype='dashed') +
    theme_light() + scale_x_continuous(breaks=seq(0, 60, 1), expand=c(0,0)) +
    ylab(TeX('$\\rho$')) + xlab('Lag')
  return(my_plot)
}

###Decorrelation function###
res_decor <- function(data, K){
  n <- nrow(data)
  index <- seq(1, n, K)
  new_ts <- data$Remainder[index]
  export <- data.frame("DATE"=data$DATE[index], "RES"=new_ts)
  return(export)
}

###plot function###
plot_decomp <- function(res, fit_df, seasons, site, variable, sp=NULL){
  #decomposition plot
  fit_df$DATE <- ymd(fit_df$DATE)
  seasons <- sort(seasons)
  nS <- length(seasons)
  end <- ncol(fit_df)
  my_range <- round(max(fit_df$Remainder), 3)
  plot_fit <- fit_df %>% pivot_longer(c(1:(3+nS),end), 
                                      names_to="key",
                                      values_to="value")
  lev_seas <- NULL
  lab_seas <- NULL
  season_levels <- for(i in 1:nS){
    lev <- paste("Seasonal", seasons[i], sep='')
    lev_seas <- c(lev_seas, lev)
    lab <- paste("Saison ", seasons[i], "j", sep='')
    lab_seas <- c(lab_seas, lab)
  }
  plot_fit$key <- factor(plot_fit$key, 
                         levels=c("Data", "Trend", 
                                  lev_seas, 
                                  "Remainder", "Predict"),
                         labels=c("Données",
                                  "Tendance",
                                  lab_seas,
                                  "Résidus", "Prédiction"))
  center <- plot_fit %>% group_by(key) %>% summarise("m"=min(value), "M"=max(value), "C"=(m+M)/2)
  scale <- data.frame("x1"=rep(ymd("2011-10-01"), length(levels(plot_fit$key))), 
                      "y1"=center$C-my_range, 
                      "y2"=center$C+my_range, 
                      "key"=levels(plot_fit$key))
  scale$key <- factor(scale$key, levels=levels(plot_fit$key))
  decomp <- plot_fit %>% ggplot(aes(DATE, value)) +
    geom_line(aes(col=key), size=.7) +
    scale_color_manual(values=c(my_palette[2:3], rep(my_palette[4], nS),
                                my_palette[6], my_palette[1])) +
    theme_light() +
    facet_wrap(key ~., scales="free_y", ncol=1, 
               strip.position = "left") +
    geom_segment(data=scale, aes(x=x1, xend=x1, y=y1, yend=y2),
                 size=3, lineend='square', col='gray60')+
    theme(legend.position="none",
          strip.text.y = element_text(size = 20),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          title=element_text(size=18)) +
    ylab(paste("log(", variable, ")", sep="")) + xlab('Temps') +
    ggtitle(paste(sp, "-", site))
  #residual plot
  res$DATE <- ymd(res$DATE)
  new_res <- res %>% ggplot(aes(x=DATE, xend=DATE, y=0, yend=RES)) + 
    geom_segment(col=my_palette[6], size=.7) +
    theme_light() + ylab("Résidus décorrélés") + xlab('Temps') +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          title=element_text(size=18))
  #residual acf
  n <- nrow(res)
  my_sd <- sd(fit_df$Remainder)
  x <- seq(min(fit_df$Remainder), max(fit_df$Remainder), 1e-4)
  y <- dnorm(x, sd=my_sd)
  df <- data.frame(x, y)
  ACF <- ggAcf(res$RES, lag.max=n) + theme_light() + 
    scale_y_continuous(limits=c(-0.2, 0.5)) + 
    ggtitle("ACF des résidus linéairement décorrélés") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          title=element_text(size=18))
  HIST <- forecast::gghistogram(fit_df$Remainder) + theme_light() + xlab("Résidus") + 
    geom_line(data=df, aes(x, y), col=my_palette[2], size=1) + 
    ylab('Compte') + ggtitle("Histogramme des résidus") +
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          title=element_text(size=18))
  #combine all plots
  res_plot <- ggarrange(NULL, new_res, ACF, NULL, HIST, NULL, nrow=6, ncol=1, heights=c(0.2, 1, 1, 0.2, 1.2, 0.2))
  my_plot <- ggarrange(decomp, NULL, res_plot, ncol=3, widths=c(1.5, 0.1, 1))
  return(my_plot)
}
  
data <- read.csv("results/TS_SYN_AB_BANYULS.csv")
plot_period(data, 'ABONDANCE', 'Banyuls', 'Abondance', 'Synechococcus', K=3)
fit <- my_mstl(data, 'ABONDANCE', 365, c(362, 183), c(11, 19), 17)
plot_cor(fit)
res <- res_decor(fit, 11)
x11()
plot_decomp(res, fit, c(362, 183), 'Banyuls', 'Abondance', 'Synechococcus')
x11()
plot(fit$Seasonal362 + fit$Seasonal183, type='l')
