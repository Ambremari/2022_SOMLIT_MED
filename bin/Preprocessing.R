#################
#SOMLIT MED - 2022
#Mathilde Couteyen Carpaye 
#Preprocessing data
################

###working directory###
setwd("C:/Users/precym-guest/Dropbox/2022_stageM2_COUTEYEN")

###packages### 
library(tidyverse)
library(naniar)
library(readr)

###data hydro###
data_hydro <- read_delim("Somlit_Extraction_hydro_20220119_101038_ea59590fd8323e08.csv",
                         delim = ";", 
                         escape_double = FALSE, 
                         trim_ws = TRUE, 
                         skip = 2)

#remove metadata
data_hydro <- data_hydro[-1,]

#missing values
data_hydro <- data_hydro %>% 
  replace_with_na_all(condition = ~.x == 999999)

#save data
write.csv(data_hydro, "2022_SOMLIT_MED/data/data_hydro.csv")

###data piconano###
data_piconano <- read_delim("Somlit_Extraction_piconano_20220103_163454_9635947c39fb31ba.csv", 
                            ";", 
                            escape_double = FALSE, 
                            trim_ws = TRUE, 
                            skip = 2)
#remove metadata
data_piconano <- data_piconano[-1,]
#missing values
data_piconano <- data_piconano %>% 
  replace_with_na_all(condition = ~.x == 999999)
data_piconano <- data_piconano %>% filter(DATE=="2012-02-14" & ID_SITE==12) %>%
  mutate(CRYSSC=NA, CRYFLR=NA)
data_piconano <- data_piconano %>% filter(DATE=="2013-01-08" & ID_SITE==11) %>%
  mutate(PICOESSC=NA, PICOEFLR=NA)

#save data
write.csv(data_hydro, "2022_SOMLIT_MED/data/data_piconano.csv")

#variable of interest
abondances <- c("TBACC", "HNABACC", "LNABACC", "CRYC", "SYNC", "PROC", "PICOEC", "NANOEC")
diff_lum <- c("TBACSSC", "HNABACSSC", "LNABACSSC", "CRYSSC", "SYNSSC", "PROSSC", "PICOESSC", "NANOESSC")

###abundance data###
ab <- data_piconano %>% dplyr:: select(abondances)
dat_ab <- cbind(data_piconano[,1:10], ab)
#wide to tidy format
dat_ab <- dat_ab %>% pivot_longer(11:18, names_to="GROUPE", values_to="ABONDANCE")
#remove NA
dat_ab <- na.omit(dat_ab)
#save data
write.csv(dat_ab, "2022_SOMLIT_MED/data/data_piconano_ab.csv")

###scatter data###
diff <- data_piconano %>% dplyr:: select(diff_lum)
dat_diff <- cbind(data_piconano[,1:10], diff)
#wide to tidy format
dat_diff <- dat_diff %>% pivot_longer(c(11:18), names_to="GROUPE", values_to="DIFFUSION")
#remove NA
dat_diff <- na.omit(dat_diff)
#save data
write.csv(dat_diff, "2022_SOMLIT_MED/data/data_piconano_diff.csv")

###data CTD###
data_CTD <-  read_delim("Somlit_Extraction_ctd_20220207_090643_e4fb36b102d70d7a.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE, 
                        skip = 2)
#remove metadata
data_CTD <- data_CTD[-1,]

###marseille###
CTD_export <- data_CTD %>% filter(ID_SITE==11)
#missing values
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999.0000)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 1e+06)
#save data
write.csv(CTD_export, "2022_SOMLIT_MED/data/data_CTD_marseille.csv")

###banyuls###
CTD_export <- data_CTD %>% filter(ID_SITE==10)
#missing values
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999.0000)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 1e+06)
#save data
write.csv(CTD_export, "2022_SOMLIT_MED/data/data_CTD_banyuls.csv")

###villefranche###
CTD_export <- data_CTD %>% filter(ID_SITE==12)
#missing values
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 999999.0000)
CTD_export <- CTD_export %>% replace_with_na_all(condition = ~.x == 1e+06)
#problematic sample
CTD_export <- CTD_export %>% filter(DATE!="2020-05-26")
#save data
write.csv(CTD_export, "2022_SOMLIT_MED/data/data_CTD_villefranche.csv")