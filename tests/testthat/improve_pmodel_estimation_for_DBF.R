#######################################################
##Aim: improve p-model performance
##specifically for the early spring and peak season
#######################################################
#As observed the main underestimation in DK-Sor(Cfb), hence in this study:
#I will first remove the data in DK-Sor, then do the calibration
#----------
library(tidyverse)
library(GenSA)
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

sites<-unique(df_recent$sitename)
for (i in 1:length(sites)) {
  df_temp<-df_recent %>%
    filter(sitename==sites[i])
  text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
  #
  df_plot<-df_temp %>%
    ggplot()+
    geom_point(aes(x=date,y=gpp))+
    geom_point(aes(x=date,y=gpp_mod),col="red")+
    annotate(geom = "text",x=text.Date,y=15,label=sites[i])
  #
  k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
  df_plot+
    geom_hline(yintercept = k,col="blue")
}

#filter the observational gpp data:
# df_recent_new<-c()
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   #filter the gpp observation data(remove the gpp that below 5 percentile and negative):
#   k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
#   df_temp$gpp[df_temp$gpp<as.numeric(k[2])& df_temp$gpp<0]<-NA
#   # df_temp %>%
#   #   ggplot()+
#   #   geom_point(aes(x=date,y=gpp))+
#   #   geom_point(aes(x=date,y=gpp_mod),col="red")+
#   #   annotate(geom = "text",x=text.Date,y=15,label=sites[i])
#   df_recent_new<-rbind(df_recent_new,df_temp)
# }
#do not filter the observation gpp in study:
df_recent_new<-df_recent

#load the data Beni sent me before:
df_old<-read.csv(file=paste0(base.path,"data/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)
#####
source(paste0(base.path,"R/","updated_R/model_hardening_byBeni_addbaseGDD.R"))
#--------------------------------------------------------------
#(2) retreive the optimized parameter for the selected sites
#--------------------------------------------------------------
# set initial value
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"k"=5)
lower=c(-50,0,0,0, 0,0)
upper=c(50,20,100,20,10,10)

# run model and compare to true values
# returns the RMSE
cost <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- model_hardening(
        .,
        par
      )

      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
    })

  df <- left_join(data, scaling_factor)

  #rmse
  # rmse <- sqrt(
  #   sum(
  #     (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  #   )/nrow(df)
  #mse:mean square error
  mse<-mean((df$gpp - df$gpp_mod * df$scaling_factor)^2,na.rm=T)
  #mae:mean absolute error:
  # mae<-sum(abs(df$gpp - df$gpp_mod * df$scaling_factor))/nrow(df)
  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)

  return(mse)
}

#--------------------------------------------------------------
#(3) optimize for each site
#--------------------------------------------------------------
#first load the PFTs information:
#load the modis data-->tidy from Beni
df_sites_modis_era <- read_rds(paste0(base.path,"data/df_sites_modis_era.csv"))
#
df_merge<-df_recent_new %>%
left_join(
    df_sites_modis_era,
    by = "sitename"
  )
##sites in Cfb-DBF
tt<-df_merge[df_merge$classid=="DBF"&df_merge$koeppen_code=="Cfb",]
unique(tt$sitename)  ##"DE-Hai" "DK-Sor" "FR-Fon" "IT-Isp"
##!!test--2022-02-25:remove the data in DK-Sor/DE-Hai and only select for DBF
df_merge<-df_merge[df_merge$sitename!="DK-Sor"&df_merge$sitename!="DE-Hai"&df_merge$sitename!="FR-Fon"
                   &df_merge$classid=="DBF",]
#main PFTs
PFTs<-unique(df_merge$classid)
# optimize for each PFT
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# par_PFTs<-c()
# for(i in 1:length(PFTs)){
#   df_sel<-df_merge %>%
#     dplyr::filter(classid==PFTs[i])
#
#   optim_par <- GenSA::GenSA(
#   par = par,
#   fn = cost,
#   data = df_sel,
#   lower = lower,
#   upper = upper,
#   control = list(max.call=5000))$par
#
#   print(i)
#   par_PFTs[[i]]<-optim_par
# }
# print("finish parameterization")
# toc()
#
# names(par_PFTs)<-PFTs
# print(par_PFTs)
# save the optimized data
# save(par_PFTs,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_forDBF.rds"))

#--------------------------------------------------------------
#(4) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
load(paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_forDBF.rds"))
#a.get the stress factor(calibration factor) for each PFT
df_final<-c()
for (i in 1:length(PFTs)) {
  df_sel<-df_merge %>%
    dplyr::filter(classid==PFTs[i])

  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening(.,par_PFTs[[i]])
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
  df_sel <- left_join(df_sel, scaling_factors)

  #merge different sites:
  df_final<-rbind(df_final,df_sel)
}

#b.make evaluation plots
#!!first need to merge the modelled gpp from different sources:
df_final$year<-lubridate::year(df_final$date)
df_merge_new<-left_join(df_final,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)
###########test for ts of temp,tmin and tmax############
#Ta
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  # group_by(sitename,doy) %>%
  summarise(gpp_obs=mean(gpp_obs_recent,na.rm=T),
          mean_Ta=mean(temp,na.rm=T),
          mean_Tmin=mean(tmin,na.rm=T),
          mean_Tmax=mean(tmax,na.rm=T),
          VPD=mean(vpd,na.rm=T),
          mean_prec=mean(prec,na.rm=T))%>%
  pivot_longer(c(mean_Ta,mean_Tmin,mean_Tmax),
               names_to = "Ta_source",values_to = "Ta") %>%
  ggplot(aes(doy,Ta,color = Ta_source))+
  geom_line()+
  facet_grid(~classid)
  # facet_grid(~sitename)
#VPD
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  summarise(VPD=mean(vpd,na.rm=T))%>%
  ggplot(aes(doy,VPD))+
  geom_line()+
  facet_grid(~classid)
#prec
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  summarise(mean_prec=mean(prec,na.rm=T))%>%
  ggplot(aes(doy,mean_prec))+
  geom_line()+
  facet_grid(~classid)
#ppfd
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  summarise(mean_ppfd=mean(ppfd,na.rm=T))%>%
  ggplot(aes(doy,mean_ppfd))+
  geom_line()+
  facet_grid(~classid)
#fapar
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  summarise(mean_fapar=mean(fapar_itpl,na.rm=T))%>%
  ggplot(aes(doy,mean_fapar))+
  geom_line()+
  facet_grid(~classid)
#gpp
df_merge_new %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid,doy) %>%
  summarise(gpp_obs=mean(gpp_obs_recent,na.rm=T),
            gpp_mod_old=mean(gpp_mod_FULL_ori,na.rm=T),
            gpp_mod_new=mean(gpp_mod_recent_ori,na.rm=T))%>%
  pivot_longer(c(gpp_obs,gpp_mod_old,gpp_mod_new),
               names_to = "gpp_source",values_to = "gpp")%>%
  ggplot(aes(doy,gpp,color=gpp_source))+
  geom_line()+
  facet_grid(~classid)


### make evaluation plots

#(1) For General plots
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(cowplot)
library(grid)

#--------------------------
#modelled and observed gpp:scatter plots
#-------------------------
plot_modobs_general<-c()
df_modobs<-c()
for(i in 1:length(PFTs)){
  #
  df_modobs_each<-df_merge_new %>%
    filter(classid==PFTs[i]) %>%
    select(sitename,date,classid,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
    mutate(gpp_obs=gpp_obs_recent,
           gpp_mod_old_ori=gpp_mod_FULL_ori,
           gpp_mod_recent_ori=gpp_mod_recent_ori,
           gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
    mutate(gpp_obs_recent=NULL,
           gpp_mod_FULL_ori=NULL)
  #
  df_modobs<-rbind(df_modobs,df_modobs_each)

  #scatter plots to compare the model and observation gpp
  gpp_modobs_comp1<-df_modobs_each %>%
    analyse_modobs2("gpp_mod_old_ori", "gpp_obs", type = "heat")
  gpp_modobs_comp2<-df_modobs_each %>%
    analyse_modobs2("gpp_mod_recent_ori", "gpp_obs", type = "heat")
  gpp_modobs_comp3<-df_modobs_each %>%
    analyse_modobs2("gpp_mod_recent_optim", "gpp_obs", type = "heat")
  # add the site-name:
  gpp_modobs_comp1$gg<-gpp_modobs_comp1$gg+
    annotate(geom="text",x=15,y=0,label=PFTs[i])
  gpp_modobs_comp2$gg<-gpp_modobs_comp2$gg+
    annotate(geom="text",x=15,y=0,label=PFTs[i])
  gpp_modobs_comp3$gg<-gpp_modobs_comp3$gg+
    annotate(geom="text",x=15,y=0,label=PFTs[i])

  #merge two plots
  evaulation_merge_plot<-plot_grid(gpp_modobs_comp1$gg,
                                   gpp_modobs_comp2$gg,gpp_modobs_comp3$gg,
                                   widths=15,heights=4,
                                   labels = "auto",ncol =3,nrow = 1,label_size = 12,align = "hv")
  # plot(evaulation_merge_plot)

  #put all the plots together:
  plot_modobs_general[[i]]<-evaulation_merge_plot
}
names(plot_modobs_general)<-PFTs

#print the plot
plot_modobs_general

#(2) For Seasonality
#a. Seasonal course for different PFTs:
#plotting:
season_plot<-df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(classid, doy) %>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="steelblue2",
                                "mod_recent_optim" = "orange", "obs" = "black"),
                     labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  annotate(geom="text",x=200,y=2,label="")+
  facet_wrap(~classid)

#print the plot
season_plot

#b. Seasonal course for each sites in different PFTs:
# For DBF:
df_modobs %>%
  filter(classid=="DBF") %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
            mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
            mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
  pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="steelblue2",
                                "mod_recent_optim" = "orange", "obs" = "black"),
                     labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  facet_wrap(~sitename)

