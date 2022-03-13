#######################################################
##Aim: improve p-model performance
##both for the early spring and peak season
#######################################################
#Try to first remove the outliers in the gpp_obs, then start to calibration
#----------
library(tidyverse)
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
#
# sites<-unique(df_recent$sitename)
# for (i in 1:length(sites)) {
#   df_temp<-df_recent %>%
#     filter(sitename==sites[i])
#   text.Date<-min(df_temp$date)+c(max(df_temp$date)-min(df_temp$date))*0.1
#   #
#   df_plot<-df_temp %>%
#     ggplot()+
#     geom_point(aes(x=date,y=gpp))+
#     geom_point(aes(x=date,y=gpp_mod),col="red")+
#     annotate(geom = "text",x=text.Date,y=15,label=sites[i])
#   #
#   k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
#   df_plot+
#     geom_hline(yintercept = k,col="blue")
# }

#filter the observational gpp data:
df_recent_new<-c()
for (i in 1:length(sites)) {
  df_temp<-df_recent %>%
    filter(sitename==sites[i])
  #filter the gpp observation data(remove the gpp that below 5 percentile and negative):
  k=quantile(df_temp$gpp,probs = c(0.01,0.05,seq(0.1,0.9,0.1)),0.95,0.99)
  df_temp$gpp[df_temp$gpp<as.numeric(k[2])& df_temp$gpp<0]<-NA
  # df_temp %>%
  #   ggplot()+
  #   geom_point(aes(x=date,y=gpp))+
  #   geom_point(aes(x=date,y=gpp_mod),col="red")+
  #   annotate(geom = "text",x=text.Date,y=15,label=sites[i])
  df_recent_new<-rbind(df_recent_new,df_temp)
}


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
sel_sites<-c("US-MMS","DK-Sor","BE-Vie","DE-Tha","NL-Loo",
             "US-UMB","US-PFa","RU-Fyo","FI-Hyy","IT-Ren")
## optimize for each site
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
par_mutisites<-c()
for(i in 1:length(sel_sites)){
  df_sel<-df_recent_new %>%
    dplyr::filter(sitename==sel_sites[i])

  optim_par <- GenSA::GenSA(
  par = par,
  fn = cost,
  data = df_sel,
  lower = lower,
  upper = upper,
  control = list(max.call=1000))$par

  print(i)
  par_mutisites[[i]]<-optim_par
}
print("finish parameterization")
toc()

names(par_mutisites)<-sel_sites
print(par_mutisites)
# save the optimized data
save(par_mutisites,file = paste0(base.path,"data/parameters_MSE_add_baseGDD/test/","optim_par_run1000_beni_mutisites.rds"))

#--------------------------------------------------------------
#(4) compare the gpp_obs, ori modelled gpp, and gpp modelled using optimated parameters
#--------------------------------------------------------------
load(paste0(base.path,"data/parameters_MSE_add_baseGDD/","optim_par_run1000_beni_mutisites.rds"))
#a.get the stress factor(calibration factor) for each site
df_final<-c()
for (i in 1:length(sel_sites)) {
  df_sel<-df_recent_new %>%
    dplyr::filter(sitename==sel_sites[i])

  scaling_factors <- df_sel %>%
    # group_by(sitename, year) %>%
    do({
      scaling_factor <- model_hardening(.,par_mutisites[[i]])
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
df_merge<-left_join(df_final,df_old,by = c("sitename", "date", "year")) %>%
  mutate(gpp_obs_recent=gpp,
         gpp_obs_old=gpp_obs,
         gpp_mod_FULL_ori=gpp_mod_FULL,
         gpp_mod_recent_ori=gpp_mod,
         gpp_mod_recent_optim=gpp_mod*scaling_factor_optim,
         gpp=NULL,
         gpp_obs=NULL,
         gpp_mod=NULL)

df_modobs<-c()
for(i in 1:length(sel_sites)){
#
df_modobs_each<-df_merge %>%
  filter(sitename==sel_sites[i]) %>%
  select(sitename,date,gpp_obs_recent,gpp_mod_FULL_ori,gpp_mod_recent_ori,gpp_mod_recent_optim) %>%
  mutate(gpp_obs=gpp_obs_recent,
         gpp_mod_old_ori=gpp_mod_FULL_ori,
         gpp_mod_recent_ori=gpp_mod_recent_ori,
         gpp_mod_recent_optim=gpp_mod_recent_optim) %>%
  mutate(gpp_obs_recent=NULL,
         gpp_mod_FULL_ori=NULL)
#
df_modobs<-rbind(df_modobs,df_modobs_each)
}


#seasonality:
sel_sites<-c("US-MMS","DK-Sor","BE-Vie","DE-Tha","NL-Loo",
             "US-UMB","US-PFa","RU-Fyo","FI-Hyy","IT-Ren")
sel_sites_type<-c("Cfa-DBF","Cfb-DBF","Cfb-MF","Cfb-ENF","Cfb-ENF",
                  "Dfb-DBF","Dfb-MF","Dfb-ENF","Dfc-ENF","Dfc-ENF")
#plotting:
#Seasonal course by climate zone:
season_plot<-df_modobs %>%
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
  annotate(geom="text",x=200,y=2,label="")+
  facet_wrap(~sitename)

