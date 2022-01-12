#refer Koen's blog for the optimization process:
#https://khufkens.com/2016/10/02/paramater-estimation-in-r-a-simple-stomatal-conductance-model-example/
# run a small example of the frost hardiness
#!!update -->using the new model that have add the Daylength(DL)
# function
library(tidyverse)
library(GenSA)
source("R/model_hardening_byBeni.R")
source("R/model_hardening_byBeni_addDL.R")

# read in data
df <- readRDS("data/model_data.rds") %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  filter(
    sitename == "US-MMS",
     date > "2008-01-01"
  ) %>%
  na.omit()

#add day length into the datasets:
library(geosphere) #calculate the day length
library(sirad)  #calculate the day of the year
library(lubridate)
#load the lat and long of the example site:
df_sites_modis_era <-read_rds(paste0("./data/df_sites_modis_era.csv"))
siteinfo<-df_sites_modis_era %>% filter(sitename=="US-MMS")
#
df$doy<-dayOfYear(df$date)
df$DL<-daylength(siteinfo$lat,df$doy)


# literature values of the
# frost hardiness function
par_old<- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1)
par_new1 <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1)
par_new2 <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1)
par_new3 <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=6)

# optimize model parameters using the
# Generalized Simulated Annealing algorithm
lower_old=c(-50,0,0,0,0)
upper_old=c(50,20,100,20,10)
#
lower_new1=c(-50,0,0,0,0)
upper_new1=c(50,20,500,20,10)
#
lower_new2=c(-50,0,0,0,0)
upper_new2=c(50,20,500,20,10)
#
lower_new3=c(-50,0,0,0,0,0)
upper_new3=c(50,20,500,20,10,24)

# run model and compare to true values
# returns the RMSE
#---
#I. for old model structure
#---
cost_old <- function(data,par_old) {
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({scaling_factor <- model_hardening(.,par_old)
      data.frame(sitename = .$sitename,
         date = .$date,
        scaling_factor = scaling_factor)
    })

  df <- left_join(data, scaling_factor)
  rmse <- sqrt(
    sum(
      (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  )/nrow(df)

  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)
  return(rmse)
}

#---
#II. for new model structure
#---
cost_new1 <- function(data,par_new1) {
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({scaling_factor <- model_hardening_new1(.,par_new1)
    data.frame(sitename = .$sitename,
               date = .$date,
               scaling_factor = scaling_factor)
    })

  df <- left_join(data, scaling_factor)
  rmse <- sqrt(
    sum(
      (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  )/nrow(df)

  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)
  return(rmse)
}

cost_new2 <- function(data,par_new2) {
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({scaling_factor <- model_hardening_new2(.,par_new2)
    data.frame(sitename = .$sitename,
               date = .$date,
               scaling_factor = scaling_factor)
    })

  df <- left_join(data, scaling_factor)
  rmse <- sqrt(
    sum(
      (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  )/nrow(df)

  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)
  return(rmse)
}

cost_new3 <- function(data,par_new3) {
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({scaling_factor <- model_hardening_new3(.,par_new3)
    data.frame(sitename = .$sitename,
               date = .$date,
               scaling_factor = scaling_factor)
    })

  df <- left_join(data, scaling_factor)
  rmse <- sqrt(
    sum(
      (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  )/nrow(df)

  # This visualizes the process,
  # comment out when running for real
  # plot(df$gpp, type = 'p',ylim=c(0,12))
  # lines(df$gpp_mod, col = "red")
  # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
  # Sys.sleep(0.1)
  return(rmse)
}

################
# optimize stuff
#-----
#I.for old model structure
#------
library(tictoc)
tic("start to parameterize")
optim_par_old <- GenSA::GenSA(
  par = par_old,
  fn = cost_old,
  data = df,
  lower = lower_old,
  upper = upper_old,
  control = list(max.call=100))$par
print("finish parameterization")
toc()

#-----
#II. for new model structure
#------
library(tictoc)
tic("start to parameterize")
optim_par_new1 <- GenSA::GenSA(
  par = par_new1,
  fn = cost_new1,
  data = df,
  lower = lower_new1,
  upper = upper_new1,
  control = list(max.call=100))$par
print("finish parameterization")
toc()

library(tictoc)
tic("start to parameterize")
optim_par_new2 <- GenSA::GenSA(
  par = par_new2,
  fn = cost_new2,
  data = df,
  lower = lower_new2,
  upper = upper_new2,
  control = list(max.call=100))$par
print("finish parameterization")
toc()

library(tictoc)
tic("start to parameterize")
optim_par_new3 <- GenSA::GenSA(
  par = par_new3,
  fn = cost_new3,
  data = df,
  lower = lower_new3,
  upper = upper_new3,
  control = list(max.call=100))$par
print("finish parameterization")
toc()

#-------
#print the parmeters
#-------
print(optim_par_old)
print(optim_par_new1)
print(optim_par_new2)
print(optim_par_new3)

###################################################
##compare the unoptimated and optimated parameters:
###################################################
#--------------------------------------------------
#a. using the optim par for the old model structure
#--------------------------------------------------
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening(.,optim_par_old)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_oldmod = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

#--------------------------------------------------
#b. using the optimilized par for the new model structure
#--------------------------------------------------
#for model structure 1
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening_new1(.,optim_par_new1)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_newmod1 = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

#for model structure 2
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening_new2(.,optim_par_new2)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_newmod2 = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

#for model structure 3
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening_new3(.,optim_par_new3)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_newmod3 = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

#----------------------
# exploratory plot
#-----------------------
p <- ggplot(df) +
  geom_point(aes(date,gpp,col="obs gpp"),col="black") +
  geom_line(aes(date,gpp_mod,col="ori model gpp")) +
  geom_line(aes(date,gpp_mod * scaling_factor_oldmod,col="optim-par old model gpp")) +
  geom_line(aes(date,gpp_mod * scaling_factor_newmod1,col="optim-par new model gpp1"))+
  geom_line(aes(date,gpp_mod * scaling_factor_newmod2,col="optim-par new model gpp2"))+
  geom_line(aes(date,gpp_mod * scaling_factor_newmod3,col="optim-par new model gpp3"))+
  labs(x = "",y = "GPP",title = "US-MMS") +
  xlim(c(as.Date("2010-01-01"),as.Date("2013-12-31")))+
  theme_minimal()

#---------------------
#model evaluaiton
#----------------------
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(cowplot)
library(grid)

#simple stats:
library(sirad)
ori_model_evas<-modeval(df$gpp_mod,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
oldmod_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_oldmod,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
newmod1_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_newmod1,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
newmod2_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_newmod2,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
newmod3_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_newmod3,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))

ori_model_evas<-t(as.data.frame(unlist(ori_model_evas)))
newmod1_par_model_evas<-t(as.data.frame(unlist(newmod1_par_model_evas)))
newmod2_par_model_evas<-t(as.data.frame(unlist(newmod2_par_model_evas)))
newmod3_par_model_evas<-t(as.data.frame(unlist(newmod3_par_model_evas)))

stats_evals<-rbind(rbind(rbind(ori_model_evas,newmod1_par_model_evas),newmod2_par_model_evas),newmod3_par_model_evas)
rownames(stats_evals)<-c("ori_model","newmod1_par_model","newmod2_par_model","newmod3_par_model")
stats_evals

#-----
#a. General evaluation:
#-----
df_modobs<-df %>%
  select(sitename,date,gpp,gpp_mod,scaling_factor_oldmod,
         scaling_factor_newmod1,scaling_factor_newmod2,scaling_factor_newmod3) %>%
  mutate(gpp_obs=gpp,
         gpp_mod_ori=gpp_mod,
         gpp_oldmod_optim=gpp_mod*scaling_factor_oldmod,
         gpp_newmod1_optim=gpp_mod*scaling_factor_newmod1,
         gpp_newmod2_optim=gpp_mod*scaling_factor_newmod2,
         gpp_newmod3_optim=gpp_mod*scaling_factor_newmod3,
         ) %>%
  mutate(gpp=NULL,
         gpp_mod=NULL)

#scatter plots to compare the model and observation gpp
gpp_modobs_comp_orimod<-df_modobs %>%
  analyse_modobs2("gpp_mod_ori", "gpp_obs", type = "hex")
gpp_modobs_comp_oldmod<-df_modobs %>%
  analyse_modobs2("gpp_oldmod_optim", "gpp_obs", type = "hex")
gpp_modobs_comp_newmod1<-df_modobs %>%
  analyse_modobs2("gpp_newmod1_optim", "gpp_obs", type = "hex")
gpp_modobs_comp_newmod2<-df_modobs %>%
  analyse_modobs2("gpp_newmod2_optim", "gpp_obs", type = "hex")
gpp_modobs_comp_newmod3<-df_modobs %>%
  analyse_modobs2("gpp_newmod3_optim", "gpp_obs", type = "hex")

#merge two plots
evaulation_merge_plot<-plot_grid(gpp_modobs_comp_orimod$gg,gpp_modobs_comp_oldmod$gg,
                                 gpp_modobs_comp_newmod1$gg,gpp_modobs_comp_newmod2$gg,
                                 gpp_modobs_comp_newmod3$gg,
                                 labels = "auto",ncol =2,label_size = 12,align = "hv")
plot(evaulation_merge_plot)

