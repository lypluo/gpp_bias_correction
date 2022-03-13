#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
#######################################################
#update in March,07:Normlize the gpp and use two parametes scaling factor( e and f)
#----------
library(tidyverse)
library(GenSA)
library(lubridate)
#-------------------------
#(1)load the data and hardening funciton
#-------------------------
base.path<-"D:/Github/gpp_bias_correction/"
#####
#load the data uploaded by Koen
df_recent <- readRDS(paste0(base.path,"data/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()
#demonstrate:
t<-df_recent %>%
  filter(sitename=="DE-Hai") %>%
  mutate(doy=lubridate::yday(date))%>%
  group_by(doy)%>%
  summarise(gpp=mean(gpp,na.rm=T))
library(phenopix)
tt<-SplineFit(t$gpp)
t_new<-data.frame(doy=t$doy,gpp=as.numeric(tt$fit$predicted))
t_new %>%
  ggplot(aes(x=doy,y=gpp))+
  geom_line()
