#######################################################
##Aim: improve p-model performance
##-->Mar,04: here I revised the dehardening functions
##both for the early spring and peak season
#######################################################
#take DK-Sor/DE-Hai/FR-Fon (has gpp underestimation in the peak season) separately to improve the model
#-->after set the iteration to 5000, the model improved substantially
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
      scaling_factor <- model_hardening_1par(
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
      scaling_factor <- model_hardening_1par(.,par_input)
      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor_optim = scaling_factor
      )
    })
  df_site_sel <- left_join(df_sel, scaling_factors)

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
#a. set initial value(one parameter:e)==> between  [spring Equinox, autumn equinox]
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"k"=5)
lower=c(-50,0,0,0,0,-10)
upper=c(50,20,200,20,2,10)
####
sel_sites<-c("DK-Sor","DE-Hai","FR-Fon")
# sel_sites<-c("FR-Fon")
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
save(par_mutisites,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run5000_beni_4several_DBFsites_rev_1parameters.rds"))

#-----------simulation demonstration-------------
plot_DK_Sor<-par_test_fun(df_recent,"DK-Sor",df_old,par_mutisites$`DK-Sor`)
plot_DE_Hai<-par_test_fun(df_recent,"DE-Hai",df_old,par_mutisites$`DE-Hai`)
plot_FR_Fon<-par_test_fun(df_recent,"FR-Fon",df_old,par_mutisites$`FR-Fon`)
#merge the plots:
library(ggpubr)
library(cowplot)
plot_grid(plot_DK_Sor,plot_DE_Hai,plot_FR_Fon,ncol = 3,align = "hv")
