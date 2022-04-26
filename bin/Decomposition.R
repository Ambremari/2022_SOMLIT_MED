#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Multi-Seasonal Trend decomposition
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN/2022_SOMLIT_MED")

###packages### 
library(tidyverse)
library(lubridate)
library(forecast)

###source code###
source("src/my_palette.R")

###MSTL function###
my_mstl <- function(data, variable, frequency, seasons, sw, tw){
  my_start <- decimal_date(ymd(data$DATE[1]))
  var <- data %>% dplyr::select(variable)
  my_ts <- msts(var, seasonal.periods=seasons,
                ts.frequency=frequency, start=my_start)
  fit <- mstl(my_ts, 
              s.window=sw, t.window=tw, iterate=2, robust=TRUE) 
  
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
  my_range <- round(max(fit_df$Remainder), 3)-0.01
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
    theme(legend.position="none") +
    ylab(paste("log(", variable, ")", sep="")) + xlab('Temps') +
    ggtitle(paste(sp, "-", site))
  #residual plot
  res$DATE <- ymd(res$DATE)
  new_res <- res %>% ggplot(aes(x=DATE, xend=DATE, y=0, yend=RES)) + 
    geom_segment(col=my_palette[6], size=.7) +
    theme_light() + ylab("Résidus décorrélés") + xlab('Temps')
  #residual acf
  n <- nrow(res)
  my_sd <- sd(res$RES)
  x <- seq(min(res$RES), max(res$RES), 1e-3)
  y <- dnorm(x, sd=my_sd)
  df <- data.frame(x, y)
  ACF <- ggAcf(res$RES, lag.max=n) + theme_light() + 
    scale_y_continuous(limits=c(-0.2, 0.5)) + ggtitle("ACF des résidus")
  HIST <- forecast::gghistogram(res$RES) + theme_light() + xlab("Résidus") + 
    geom_line(data=df, aes(x, y), col=my_palette[2], size=1) + 
    ylab('Compte') + ggtitle("Histogramme des résidus")
  #combine all plots
  res_plot <- ggarrange(NULL, new_res, ACF, HIST, NULL, nrow=5, ncol=1, heights=c(0.5, 1, 1, 1.2, 0.5))
  my_plot <- ggarrange(decomp, NULL, res_plot, ncol=3, widths=c(1.5, 0.1, 1))
  return(my_plot)
}
  

