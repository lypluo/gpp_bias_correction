#----------------------------------------------------
# run a small example of the frost hardiness function
# compare the method adopt from Mekela et al., 2008 and Beni's method
#Aim: to confirm and test if adding fT could improve the model efficiency
#-->results: I found it takes long time to parameterization for fT_fL-->hence,
#to plan to only use fT...
#----------------------------------------------------
library(tidyverse)
library(GenSA)
library(lubridate)
# source("R/model_hardening_byBeni.R")
source("R/updated_R/new_formulated/model_fT.R")
source("R/updated_R/new_formulated/model_fT_rev.R")
source("R/updated_R/new_formulated/model_fT_rev_fL.R")
#-------------
#(1)read in data
#-------------
df <- readRDS("data/model_data.rds") %>%
  mutate(year = format(date, "%Y"),
         doy =yday(date)
         # ppfd=ppfd*1000000
         ) %>%
  filter(
    sitename == "US-Ha1",
    date > "2010-01-01"
  )

#----------------------------
#(2)retreive the optimized parameter
#----------------------------
# set initial value
# par_Beni <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1)
par_YP <- c("tau"=5,"X0"=-10,"Smax"=5)
par_YP_rev <- c("tau"=5,"X0"=-10,"Smax"=5,"k"=1)
par_YP_fTfL<- c("tau"=5,"X0"=-10,"Smax"=5,"lamda"=0,"k"=1)
# optimize model parameters using the Generalized Simulated Annealing algorithm
#refer Koen's blog for the optimization process:
#https://khufkens.com/2016/10/02/paramater-estimation-in-r-a-simple-stomatal-conductance-model-example/
#set parameter variation range-->refer Tian et al., 2021
# lower_Beni=c(-50,0,0,0,0)
# upper_Beni=c(50,20,100,20,10)
lower_YP=c(1,-10,5)
upper_YP=c(25,10,25)
lower_YP_rev=c(1,-10,5,0)
upper_YP_rev=c(25,10,25,2)
lower_YP_fTfL<-c(1,-10,5,0,0)
upper_YP_fTfL<-c(25,10,25,0.5,2)


# run model and compare to true values
# returns the RMSE
#for Beni:
# cost_Beni <- function(
#   data,
#   par
# ) {
#
#   scaling_factor <- data %>%
#     # group_by(sitename) %>%
#     do({
#       scaling_factor <- model_hardening(
#         .,
#         par
#       )
#
#       data.frame(
#         sitename = .$sitename,
#         date = .$date,
#         scaling_factor = scaling_factor
#       )
#     })
#
#   df <- left_join(data, scaling_factor)
#   #remove the NA lines in df:
#   df<-df[!is.na(df$gpp),]
#
#   rmse <- sqrt(
#     sum(
#       (df$gpp - df$gpp_mod * df$scaling_factor)^2)
#   )/nrow(df)
#
#   # This visualizes the process,
#   # comment out when running for real
#   # plot(df$gpp, type = 'p',ylim=c(0,12))
#   # lines(df$gpp_mod, col = "red")
#   # lines(df$gpp_mod * df$scaling_factor, col = "blue",cex=1.2)
#   # Sys.sleep(0.1)
#
#   return(rmse)
# }
#for YP:
###For 3 pars:
cost_YP <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- f_Ts(
        .,
        par,
        plot = FALSE
      )

      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
    })

  df <- left_join(data, scaling_factor)
  #remove the NA lines in df:
  df<-df[!is.na(df$gpp),]

  # rmse <- sqrt(
  #   sum(
  #     (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  # )/nrow(df)
  #mse:mean square error
  mse<-mean((df$gpp - df$gpp_mod * df$scaling_factor)^2,na.rm=T)
  return(mse)
}

###For 4 pars:
cost_YP_rev <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- f_Ts_rev(
        .,
        par,
        plot = FALSE
      )

      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
    })

  df <- left_join(data, scaling_factor)
  #remove the NA lines in df:
  df<-df[!is.na(df$gpp),]

  # rmse <- sqrt(
  #   sum(
  #     (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  # )/nrow(df)
  #mse:mean square error
  mse<-mean((df$gpp - df$gpp_mod * df$scaling_factor)^2,na.rm=T)
  return(mse)
}

###For pars for fTfL-->the parameterization took too much time-->
#at the end, do not plan to add fL..
cost_YP_fTfL <- function(
  data,
  par
) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- f_Ts_rev_fL(
        .,
        par,
        plot = FALSE
      )

      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
    })

  df <- left_join(data, scaling_factor)
  #remove the NA lines in df:
  df<-df[!is.na(df$gpp),]

  # rmse <- sqrt(
  #   sum(
  #     (df$gpp - df$gpp_mod * df$scaling_factor)^2)
  # )/nrow(df)
  #mse:mean square error
  mse<-mean((df$gpp - df$gpp_mod * df$scaling_factor)^2,na.rm=T)
  return(mse)
}

#-----------------
# optimize stuff
#-----------------
#for Beni's function
# library(tictoc)#-->record the parameterization time
# tic("start to parameterize")
# optim_par_Beni <- GenSA::GenSA(
#   par = par_Beni,
#   fn = cost_Beni,
#   data = df,
#   lower = lower_Beni,
#   upper = upper_Beni,
#   control = list(max.call=10))$par
# print("finish parameterization")
# toc()

#for YP's function
#for 3 pars
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
optim_par_YP_3pars <- GenSA::GenSA(
  par = par_YP,
  fn = cost_YP,
  data = df,
  lower = lower_YP,
  upper = upper_YP,
  control = list(max.call=1000))$par
print("finish parameterization")
toc()

#for 4 pars
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
optim_par_YP_4pars <- GenSA::GenSA(
  par = par_YP_rev,
  fn = cost_YP_rev,
  data = df,
  lower = lower_YP_rev,
  upper = upper_YP_rev,
  control = list(max.call=1000))$par
print("finish parameterization")
toc()

#for pars for fTfL:5 pars
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
optim_par_YP_5pars <- GenSA::GenSA(
  par = par_YP_fTfL,
  fn = cost_YP_fTfL,
  data = df,
  lower = lower_YP_fTfL,
  upper = upper_YP_fTfL,
  control = list(max.call=2))$par
print("finish parameterization")
toc()

#--------------
#(3)get the stress factor using the optimized paramters:
#--------------
#for Beni:
# scaling_factors_Beni <- df %>%
#   group_by(sitename, year) %>%
#   do({
#     scaling_factor <- model_hardening(
#       .,
#       optim_par_Beni,
#       plot = T
#     )
#
#     data.frame(
#       sitename = .$sitename,
#       date = .$date,
#       scaling_factor = scaling_factor
#     )
#   })
#
# df_Beni <- left_join(df, scaling_factors_Beni)

###
#for YP:
#3pars
scaling_factors_YP_3pars <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- f_Ts(
      .,
      optim_par_YP_3pars,
      plot = T
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor = scaling_factor
    )
  })

df_YP_3pars <- left_join(df, scaling_factors_YP_3pars)
#4pars:
scaling_factors_YP_4pars <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- f_Ts_rev(
      .,
      optim_par_YP_4pars,
      plot = T
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor = scaling_factor
    )
  })

df_YP_4pars <- left_join(df, scaling_factors_YP_4pars)

#5pars:
scaling_factors_YP_5pars <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- f_Ts_fTfL(
      .,
      optim_par_YP_5pars,
      plot = T
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor = scaling_factor
    )
  })

df_YP_5pars <- left_join(df, scaling_factors_YP_5pars)

#-----------------
##merge the data
#-----------------
# df_Beni<-df_Beni %>%
#   mutate(GPP_Beni=gpp_mod*scaling_factor,
#          scaling_factor=NULL)
df_YP_3pars<-df_YP_3pars %>%
  mutate(gpp_YP_3pars=gpp_mod*scaling_factor,
         scaling_factor=NULL)
df_YP_4pars<-df_YP_4pars %>%
  mutate(gpp_YP_4pars=gpp_mod*scaling_factor,
         scaling_factor=NULL)
df_YP_5pars<-df_YP_4pars %>%
  mutate(gpp_YP_5pars=gpp_mod*scaling_factor,
         scaling_factor=NULL)

# df_merge<-cbind(df_Beni,df_YP[,"GPP_YP"])
df_merge<-cbind(df_YP_3pars,
                df_YP_4pars[,"gpp_YP_4pars"],
                df_YP_5pars[,"gpp_YP_5pars"])
#-------------------
#plotting
#-------------------
###################
# exploratory plot
# p <- ggplot(df_merge) +
#   geom_line(
#     aes(date,GPP_Beni),colour = "red") +
#   geom_line(
#     aes(date,GPP_YP),colour = "blue") +
#   geom_point(
#     aes(date,gpp)) +
#   labs(
#     x = "",
#     y = "GPP",
#     title = "US-NR1"
#   ) +
#   theme_minimal()

p <- ggplot(df_merge) +
  geom_line(
    aes(date,gpp_mod),colour = "red") +
  # geom_line(
  #   aes(date,gpp_YP_3pars),colour = "blue") +
  # geom_line(
  #     aes(date,gpp_YP_4pars),colour = "orange") +
  # geom_line(
  #     aes(date,gpp_YP_5pars),colour = "grey") +
  geom_point(
    aes(date,gpp)) +
  labs(
    x = "",
    y = "GPP",
    title = "US-NR1"
  ) +
  theme_minimal()

print(p)

# ggsave("niwot.png", width = 7, height = 4)
########################
#Evaluation: scattering plot using functions in Beni's package:
#----------------
library(rbeni)
library(ggpubr)
library(cowplot)
library(grid)
p_ori_mod<-df_merge %>%
analyse_modobs2("gpp_mod", "gpp", type = "heat")
p_rev_mod1<-df_merge %>%
  analyse_modobs2("gpp_YP_3pars", "gpp", type = "heat")
p_rev_mod2<-df_merge %>%
  analyse_modobs2("gpp_YP_4pars", "gpp", type = "heat")
#
plot_grid(p_ori_mod,p_rev_mod1,p_rev_mod2,nrow = 1)

