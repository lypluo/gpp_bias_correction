##---------------------------------------
#Aim: To compare the parameters among different groups(e.g.PFTs) after calibrating
#parameters for each site
##---------------------------------------
library(dplyr)
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(tidyverse)
#---------------------------
#(1)load the calibrated parameters for each site
#---------------------------
base.path<-"D:/Github/gpp_bias_correction/"
load(paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_eachsite_updated.rds"))
#merge the parameters:
merge_pars<-c()
sites<-names(par_mutisites)
for(i in 1:length(par_mutisites)){
  temp<-t(as.data.frame(par_mutisites[i]))
  merge_pars<-rbind(merge_pars,temp)
}
pars_final<-as.data.frame(t(merge_pars))
names(pars_final)<-sites

#----------------------------
#(2)load original data(meteos, gpp..)
#----------------------------
#load the data uploaded by Koen
df_recent <- readRDS(paste0(base.path,"data/model_data.rds")) %>%
  mutate(
    year = format(date, "%Y")
  )
#load the PFTs information:
#load the modis data-->tidy from Beni
#read_rds from tidyverse
df_sites_modis_era <- read_rds(paste0(base.path,"data/df_sites_modis_era.csv"))
#----
#merge the data
#-----
df_merge<-left_join(df_recent,df_sites_modis_era,by="sitename")
df_merge$year<-as.numeric(df_merge$year)
#load the data Beni sent me before:
df_old<-read.csv(file=paste0(base.path,"data/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date))
#----
#merge data:
#----
df_merge_new<-left_join(df_merge,df_old,by=c("sitename", "date", "year"))

#-----------
#summarize the meteos for each site
#-----------
df_sum_yearly_1<-df_merge_new %>%
  group_by(sitename,year) %>%
  summarise(temp=mean(temp),
            prec=sum(prec),
            vpd=mean(vpd),
            ppdf=mean(ppfd),
            elv=mean(elv),
            tmin=mean(tmin),
            tmax=mean(tmax),
            fapar_itpl=mean(fapar_itpl),
            fapar_spl=mean(fapar_spl)
            )
df_sum_yearly_2<-df_merge_new %>%
  group_by(sitename,year) %>%
  summarise(lon=unique(lon),
            lat=unique(lat),
            classid=unique(classid),
            koeppen_code=unique(koeppen_code))
df_sum_yearly<-left_join(df_sum_yearly_1,df_sum_yearly_2)

#---
#summary site-year for site
#---
df_sum_1<-df_sum_yearly %>%
  group_by(sitename) %>%
  summarise_at(vars(temp:fapar_spl),mean,na.rm=T)
df_sum_2<-df_sum_yearly %>%
  group_by(sitename) %>%
  summarise(lon=unique(lon),
            lat=unique(lat),
            classid=unique(classid),
            koeppen_code=unique(koeppen_code))
df_sum<-left_join(df_sum_1,df_sum_2)

##-----------------------
#(3) compare the parameter difference in differnt group
##----------------------
df_sum$Clim.PFTs<-paste0(df_sum$koeppen_code,"-",df_sum$classid)
#only target the sites we used for the analysis:
##---------------------
#A.load the event_length data-->the sites we were used
#---------------------
load.path<-"D:/data/photocold_project/event_length/Using_sites_in_Fluxnet2015/"
load(paste0(load.path,"df_events_length.RDA"))
#
used_sites<-unique(df_events_all$sitename)

#-----select the data for thoses used sites----
df_final<-df_sum %>%
  filter(sitename %in% used_sites)

#--------------------------
#(4)plotting:
#--------------------------
#first merge the parameters with meteos:
pars_final<-as.data.frame(t(pars_final))
pars_final$sitename<-rownames(pars_final)
#merge:
df_final_new<-left_join(df_final,pars_final,by="sitename")
#---
#check the variables distribution and boxplots
#---
vars.names<-c("a","b","c","d","e","f","k")
for (i in 1:length(vars.names)) {
  hist(as.numeric(unlist(df_final_new[,vars.names[i]])),xlab = vars.names[i])
}
#boxplot:
data_sel<-df_final_new %>%
  select(sitename,classid,a:k)%>%
  pivot_longer(c(a:k),names_to = "parameter",values_to = "parameter_value")
ggplot(data=data_sel,aes(x=parameter,y=parameter_value))+
  geom_boxplot()+
  geom_jitter(aes(col=classid),width = 0.1)+
  facet_wrap(~parameter,scales = "free")


#----
#check the paraters difference among different groups
#----
library(ggpubr)
library(cowplot)
check_groups<-function(df,par_name){
  # df<-df_final_new
  # par_name<-"a"

  df_t<-df %>%
    select(sitename,classid,koeppen_code,Clim.PFTs,par_name)
  names(df_t)<-c("sitename","classid","koeppen_code","Clim.PFTs","par")
  #for different PFTs
  p_PFTs<-ggplot(data=df_t,aes(x=par,color=classid,fill=classid))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    geom_density(alpha=.2)+
    xlab(par_name)
  #for different Clim.
  p_Clim<-ggplot(data=df,aes(x=b,color=koeppen_code,fill=koeppen_code))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    geom_density(alpha=.2)+
    xlab(par_name)
  #for different Clim.-PFTs
  p_Clim.PFTs<-ggplot(data=df,aes(x=b,color=Clim.PFTs,fill=Clim.PFTs))+
    geom_histogram(aes(y=..density..,),position = "identity",binwidth = 1,alpha=0.5)+
    xlab(par_name)
  # geom_density(alpha=.2)
  #
  p_merge<-plot_grid(p_PFTs,p_Clim,p_Clim.PFTs)
  return(p_merge)
}
##
check_groups(df_final_new,"a")
check_groups(df_final_new,"b")
check_groups(df_final_new,"c")
check_groups(df_final_new,"d")
check_groups(df_final_new,"e")
check_groups(df_final_new,"f")
check_groups(df_final_new,"k")

#----
#check environmental drivers relationship between parameters
#----
check_relation<-function(df,par_name){
  # df<-df_final_new
  # par_name<-"a"

  df_t<-df %>%
    select(sitename:fapar_spl,classid,koeppen_code,Clim.PFTs,par_name)
  names(df_t)<-c("sitename","temp","prec","vpd","ppdf","elv",
                 "tmin","tmax","fapar_itpl","fapar_spl",
                 "classid","koeppen_code","Clim.PFTs","par")
  #
  p_ta<-ggplot(data=df_t,aes(x=temp,y=par,color=classid))+
    geom_point()+
    xlab("ta")+
    ylab(par_name)
  #
  p_prec<-ggplot(data=df_t,aes(x=prec,y=par,color=classid))+
    geom_point()+
    xlab("prec")+
    ylab(par_name)
  p_vpd<-ggplot(data=df_t,aes(x=vpd,y=par,color=classid))+
    geom_point()+
    xlab("vpd")+
    ylab(par_name)
  p_ppfd<-ggplot(data=df_t,aes(x=ppdf,y=par,color=classid))+
    geom_point()+
    xlab("ppfd")+
    ylab(par_name)
  p_tmin<-ggplot(data=df_t,aes(x=tmin,y=par,color=classid))+
    geom_point()+
    xlab("tmin")+
    ylab(par_name)
  p_tmax<-ggplot(data=df_t,aes(x=tmax,y=par,color=classid))+
    geom_point()+
    xlab("tmax")+
    ylab(par_name)
  p_fapar<-ggplot(data=df_t,aes(x=fapar_itpl,y=par,color=classid))+
    geom_point()+
    xlab("fapar_itpl")+
    ylab(par_name)
    #
  p_merge<-plot_grid(p_ta,p_prec,p_vpd,
                     p_ppfd,p_tmin,p_fapar,nrow = 2,align = "hv")
  return(p_merge)
}

#
check_relation(df_final_new,"a")
check_relation(df_final_new,"b")
check_relation(df_final_new,"c")
check_relation(df_final_new,"d")
check_relation(df_final_new,"e")
check_relation(df_final_new,"f")
check_relation(df_final_new,"k")

