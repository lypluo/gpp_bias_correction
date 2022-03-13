#######################################################
##Aim: improve p-model performance
##-->Mar,04: here I revised the dehardening functions
##both for the early spring and peak season
#March,07,normalization the gpp
#######################################################
#----------
library(tidyverse)
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
#add doy
df_recent$doy<-lubridate::yday(df_recent$date)
#
# sites<-unique(df_recent$sitename)
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
#   df_temp %>%
#     ggplot()+
#     geom_point(aes(x=date,y=gpp))+
#     geom_point(aes(x=date,y=gpp_mod),col="red")+
#     annotate(geom = "text",x=text.Date,y=15,label=sites[i])
# }

#first load the PFTs information:
#load the modis data-->tidy from Beni
df_sites_modis_era <- read_rds(paste0(base.path,"data/df_sites_modis_era.csv"))
#
df_recent<-df_recent %>%
  left_join(
    df_sites_modis_era,
    by = "sitename"
  )
#day of the year
library(sirad)  #calculate the day of the year
df_recent$doy<-dayOfYear(df_recent$date)
#identify the Cfa-DBF sites:
#main Clim-PFTs
df_recent$Clim_PFTs<-paste0(df_recent$koeppen_code,"-",df_recent$classid)
unique(df_recent[df_recent$Clim_PFTs=="Cfa-DBF",]$sitename)

#------------------------------------------
#(3)normalized the GPP-->for each site,
#normalized the "gpp" and "gpp_mod" through their 90 percentiles
#------------------------------------------
#I should use the same value to normlize the gpp and gpp_mod:
gpp_P95<-df_recent %>%
  group_by(sitename) %>%
  summarise(gpp_norm_p95=quantile(c(gpp,gpp_mod),0.95,na.rm=T))
#
df_recent<-left_join(df_recent,gpp_P95,by="sitename")
df_recent<-df_recent %>%
  mutate(gpp=gpp/gpp_norm_p95,gpp_mod=gpp_mod/gpp_norm_p95)


#load the data Beni sent me before:
df_old<-read.csv(file=paste0(base.path,"data/","ddf_fluxnet2015_pmodel_with_forcings_stocker19gmd.csv"))
df_old<-df_old %>%
  mutate(date=lubridate::mdy(date),
         year=lubridate::year(date)) %>%
  na.omit(gpp_obs)
#####
source(paste0(base.path,"R/","updated_R/model_hardening_byBeni_addbaseGDD_rev.R"))
# source(paste0(base.path,"R/","updated_R/model_hardening_byBeni_addbaseGDD.R"))
#--------------------------------------------------------------
#(2) retreive the optimized parameter for the selected sites
#--------------------------------------------------------------
# set initial value
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=1,"k"=5)
lower=c(-50,0,0,0,0,0,0)
upper=c(50,20,100,20,2,2,10)

# run model and compare to true values
# returns the RMSE
cost <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- model_hardening_2par(
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

###-----------------
#(3)Further parameter test for the site DK-Sor:
###----------------
#first to check the optimilized parameters:
# par_optimized<-par_mutisites$`DK-Sor`
#
par_test_fun<-function(df_recent,site,df_old,par_input){
  # df_recent<-df_recent
  # site<-"DK-Sor"
  # df_old<-df_old
  # par_input<-par_optimized


  ##--------selecting the data-----------
  df_site_sel<-c()
  df_sel<-df_recent %>%
    dplyr::filter(sitename==site)

  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening_2par(.,par_input)
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
  df_site_sel <- left_join(df_sel, scaling_factors)
  #-----------------------------
  #update:March-07
  #need to back-convert the normalized gpp to gpp
  #-----------------------------
  df_site_sel<-df_site_sel %>%
    mutate(gpp=gpp*gpp_norm_p95,
           gpp_mod=gpp_mod*gpp_norm_p95)

  df_site_sel$year<-lubridate::year(df_site_sel$date)
  df_merge<-left_join(df_site_sel,df_old,by = c("sitename", "date", "year")) %>%
    mutate(gpp_obs_recent=gpp,
           gpp_obs_old=gpp_obs,
           gpp_mod_FULL_ori=gpp_mod_FULL,
           gpp_mod_recent_ori=gpp_mod,
           gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
           gpp=NULL,
           gpp_obs=NULL,
           gpp_mod=NULL)
  df_modobs<-c()
  ##------mege the data-----------------
  df_modobs_each<-df_merge %>%
    filter(sitename==site) %>%
    select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
    mutate(gpp_obs=gpp_obs_recent,
           gpp_mod_old_ori=gpp_mod_FULL_ori,
           gpp_mod_recent_ori=gpp_mod_recent_ori,
           gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
    mutate(gpp_obs_recent=NULL,
           gpp_mod_FULL_ori=NULL)
  ##plotting---------------------------
  season_plot<-df_modobs_each %>%
    mutate(doy = lubridate::yday(date)) %>%
    group_by(sitename, doy) %>%
    summarise(obs = mean(gpp_obs, na.rm = TRUE),
              mod_old_ori=mean(gpp_mod_old_ori, na.rm = TRUE),
              mod_recent_ori=mean(gpp_mod_recent_ori, na.rm = TRUE),
              mod_recent_optim=mean(gpp_mod_recent_optim,na.rm = TRUE)) %>%
    pivot_longer(c(obs,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
    # pivot_longer(c(obs,mod_old_ori,mod_recent_ori,mod_recent_optim), names_to = "Source", values_to = "gpp") %>%
    ggplot(aes(doy, gpp, color = Source)) +
    geom_line() +
    scale_color_manual(values = c("mod_old_ori" = "red","mod_recent_ori"="steelblue2",
                                  "mod_recent_optim" = "orange", "obs" = "black"),
                       labels = c("Old P-model","Recent Ori P-model", "Recent Optim P-model","Obs.")) +
    labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
         x = "Day of year") +
    annotate(geom="text",x=200,y=2,label="")+
    facet_wrap(~sitename)

  #
  print(season_plot)
  return(season_plot)
}

##---------------------------
##(4)re-set parameter range and calibrate":
##---------------------------
#b. set initial value(two parameters:e and f)
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=1,"k"=5)
lower=c(-50,0,0,0,0,0,-10)
upper=c(50,20,200,20,2,2,10)
####
#selected the Cfa-DBF sites:
sel_sites<-c("IT-Col","IT-PT1","US-MMS")
# sel_sites<-c("FR-Fon",
##optimize for each site
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
par_mutisites<-c()
for(i in 1:length(sel_sites)){
  df_sel<-df_recent %>%
    dplyr::filter(sitename==sel_sites[i])

  optim_par <- GenSA::GenSA(
  par = par,
  fn = cost,
  data = df_sel,
  lower = lower,
  upper = upper,
  control = list(max.call=5000))$par

  print(paste0("===============",i,"================"))
  par_mutisites[[i]]<-optim_par
}
print("finish parameterization")
toc()
#
names(par_mutisites)<-sel_sites
par_mutisites<-par_mutisites
# save the optimized data
# save(par_mutisites,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_4several_Cfa_DBF_rev_2parameters.rds"))
load(paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_4several_Cfa_DBF_rev_2parameters.rds"))

#-----------simulation demonstration-------------
plot_IT_Col1<-par_test_fun(df_recent,"IT-Col",df_old,par_mutisites$`IT-Col`)
plot_IT_PT11<-par_test_fun(df_recent,"IT-PT1",df_old,par_mutisites$`IT-PT1`)
plot_US_MMS1<-par_test_fun(df_recent,"US-MMS",df_old,par_mutisites$`US-MMS`)
#merge the plots:
library(ggpubr)
library(cowplot)
plot_grid(plot_IT_Col1,plot_IT_PT11,plot_US_MMS1,ncol = 3,align = "hv")

#---------------
#additional:
#---------------
#using the same parameters for Cfa-DBF
load(paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_Clim_andPFTs_update.rds"))
plot_IT_Col2<-par_test_fun(df_recent,"IT-Col",df_old,par_Clim_PFTs$`Cfa-DBF`)
plot_IT_PT12<-par_test_fun(df_recent,"IT-PT1",df_old,par_Clim_PFTs$`Cfa-DBF`)
plot_US_MMS2<-par_test_fun(df_recent,"US-MMS",df_old,par_Clim_PFTs$`Cfa-DBF`)
#merge the plots:
plot_grid(plot_IT_Col1,plot_IT_Col2,ncol = 2,align = "hv")
plot_grid(plot_US_MMS1,plot_US_MMS2,ncol = 2,align = "hv")
