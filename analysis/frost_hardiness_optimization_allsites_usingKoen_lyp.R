#-----------------------------------------------------
#get the general optimized parameter for all sites
#Refer Koen's optimization script==>update this by YP on Dec,3, 2021
#-----------------------------------------------------
library(tidyverse)
library(GenSA)
source("R/frost_hardiness.R")

#--------------------------
#(1)read data==>updated by YP: using all the site-years data
#--------------------------
df <- readRDS("data/model_data.rds") %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()
#----------------------------
#(2)retreive the optimized parameter
#----------------------------
# literature values of the frost hardiness function (initial value)
par <- c(
  a = 1.5,
  b = -21.5,
  T8 = 11.3,
  t = 5,
  base_t = -25
)

# optimize model parameters using the Generalized Simulated Annealing algorithm
#refer Koen's blog for the optimization process:
#https://khufkens.com/2016/10/02/paramater-estimation-in-r-a-simple-stomatal-conductance-model-example/
#set parameter variation range:
lower=c(0,-100,-50,1, -100)
upper=c(100,0,100,90, 0)

# run model and compare to true values
# returns the RMSE
cost <- function(
  data,
  par
  ) {

  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({
      scaling_factor <- frost_hardiness(
        #.$temp,                  #update by YP: using the Tmin according to literature
        .$tmin,
        par
      )

      data.frame(
        sitename = .$sitename,
        date = .$date,
        scaling_factor = scaling_factor
      )
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

# optimize stuff
library(tictoc)#-->record the parameterization time
tic("start to parameterize")
optim_par <- GenSA::GenSA(
  par = par,
  fn = cost,
  data = df,
  lower = lower,
  upper = upper,
  control = list(max.call=100))$par
print("finish parameterization")
toc()
#save the optimized data
save(optim_par,file = paste0("./data/parameters/","optim_par_run100_koen.rds"))
# load(paste0("./data/parameters/","optim_par_run100_koen.rds"))

#----------------------------
##(3)compare the results using defaulted and optimated parameters-->site by site
#----------------------------
#a. using the default par
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- frost_hardiness(.$tmin,par)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_default = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)
#b. using the optimilized par
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- frost_hardiness(.$tmin,optim_par)
    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_optim = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

# exploratory plot
site_names<-unique(df$sitename)
p_all<-c()
for(i in 1:length(site_names)){
  site_run<-site_names[i]
  df_run<-df %>% filter(sitename==site_run)
  label_x<-min(df_run$date)+91
  label_y<-round(max(df_run$gpp)+1,0)
  p <- ggplot(df_run) +
    geom_point(aes(date,gpp,col="obs"),size=2) +
    geom_line(aes(date,gpp_mod,col="ori-mod")) +
    geom_line(aes(date,gpp_mod * scaling_factor_default,col="default-par model gpp")) +
    geom_line(aes(date,gpp_mod * scaling_factor_optim,col="optim-par model gpp"))+
    labs(x = "",y = "GPP") +
    scale_color_manual(values = c("ori-mod"=adjustcolor("tomato",0.8),"default-par model gpp"=adjustcolor("steelblue2",0.8),
                                  "optim-par model gpp"=adjustcolor("green3",0.8),"obs"=adjustcolor("grey",0.8)))+
    annotate(geom = "text",x=label_x,y=label_y,label=site_run)+
    # facet_grid(sitename ~ .) +
    theme_minimal()
  p_all[[i]]<-p
}
names(p_all)<-site_names

#---------------------------
##(4)model evaluation:
#---------------------------
#simple stats:
library(sirad)
ori_model_evas<-modeval(df$gpp_mod,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
default_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_default,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))
optim_par_model_evas<-modeval(df$gpp_mod*df$scaling_factor_optim,df$gpp,stat=c("N","pearson","MAE","RMSE","R2","slope","intercept","EF"))

ori_model_evas<-t(as.data.frame(unlist(ori_model_evas)))
default_par_model_evas<-t(as.data.frame(unlist(default_par_model_evas)))
optim_par_model_evas<-t(as.data.frame(unlist(optim_par_model_evas)))

stats_evals<-rbind(rbind(ori_model_evas,default_par_model_evas),optim_par_model_evas)
rownames(stats_evals)<-c("ori_model","default_par_model","optim_par_model")

#evaluation using Beni's function:
devtools::load_all("D:/Github/rbeni/")
library(rbeni) #-->make the evaluation plot
library(cowplot)
library(grid)

##(2)make the evaluation plot
#-----
#a. General evaluation:
#-----
df_modobs<-df %>%
  select(sitename,date,gpp,gpp_mod,scaling_factor_default,scaling_factor_optim) %>%
  mutate(gpp_obs=gpp,
         gpp_mod_ori=gpp_mod,
         gpp_mod_default_par=gpp_mod*scaling_factor_default,
         gpp_mod_optim_par=gpp_mod*scaling_factor_optim) %>%
  mutate(gpp=NULL,
         gpp_mod=NULL)

#scatter plots to compare the model and observation gpp
gpp_modobs_comp1<-df_modobs %>%
  analyse_modobs2("gpp_mod_ori", "gpp_obs", type = "hex")
gpp_modobs_comp2<-df_modobs %>%
  analyse_modobs2("gpp_mod_default_par", "gpp_obs", type = "hex")
gpp_modobs_comp3<-df_modobs %>%
  analyse_modobs2("gpp_mod_optim_par", "gpp_obs", type = "hex")

#merge two plots
evaulation_merge_plot<-plot_grid(gpp_modobs_comp1$gg,
    gpp_modobs_comp2$gg,gpp_modobs_comp3$gg,
labels = "auto",ncol =3,label_size = 12,align = "hv")
plot(evaulation_merge_plot)

### Seasonality
#-------------
#b.By climate zone
#------------
#load the modis data-->tidy from Beni
df_sites_modis_era <- read_rds(paste0("./data/df_sites_modis_era.csv"))
#now take the df_modobs(ori,default, and optimized model gpp as the comparison) as the example
df_season_climate <- df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_ori = mean(gpp_mod_ori, na.rm = TRUE),
            mod_default = mean(gpp_mod_default_par, na.rm = TRUE),
            mod_optim = mean(gpp_mod_optim_par, na.rm = TRUE)) %>%
  left_join(
    df_sites_modis_era,
    by = "sitename"
  ) %>%
  mutate(northsouth = ifelse(lat>0, "North", "South")) %>%
  dplyr::filter(koeppen_code != "-") %>%
  mutate(kg_code_northsouth = paste(koeppen_code, northsouth)) %>%
  group_by(kg_code_northsouth, doy) %>%
  summarise(obs = mean(obs, na.rm = TRUE),
            mod_ori = mean(mod_ori, na.rm = TRUE),
            mod_default = mean(mod_default, na.rm =T),
            mod_optim = mean(mod_optim, na.rm = TRUE)
            )

#plotting:
#Seasonal course by climate zone:
df_season_climate %>%
  pivot_longer(c(obs, mod_ori,mod_default,mod_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_ori" = "red","mod_default" = "green3",
                                "mod_optim" = "orange", "obs" = "black"),
        labels = c("Ori P-model","Default-par P-model","Optim-par P-model", "obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  facet_wrap(~kg_code_northsouth)

ggsave("./fig/meanseasonalcycle_by_climate_cal_fromKoen.pdf", width = 9, height = 6)

#-------------
#b.By each site
#------------
#same as by the climate
df_season_site <- df_modobs %>%
  mutate(doy = lubridate::yday(date)) %>%
  group_by(sitename, doy) %>%
  summarise(obs = mean(gpp_obs, na.rm = TRUE),
            mod_ori = mean(gpp_mod_ori, na.rm = TRUE),
            mod_default = mean(gpp_mod_default_par, na.rm = TRUE),
            mod_optim = mean(gpp_mod_optim_par, na.rm = TRUE)) %>%
  group_by(sitename, doy) %>%
  summarise(obs = mean(obs, na.rm = TRUE),
            mod_ori = mean(mod_ori, na.rm = TRUE),
            mod_default = mean(mod_default, na.rm =T),
            mod_optim = mean(mod_optim, na.rm = TRUE)
  )

#plotting:
#Seasonal course by climate zone:
df_season_site %>%
  pivot_longer(c(obs, mod_ori,mod_default,mod_optim), names_to = "Source", values_to = "gpp") %>%
  ggplot(aes(doy, gpp, color = Source)) +
  geom_line() +
  scale_color_manual(values = c("mod_ori" = "red","mod_default" = "green3",
                                "mod_optim" = "orange", "obs" = "black"),
                     labels = c("Ori P-model","Default-par P-model","Optim-par P-model", "obs.")) +
  labs(y = expression( paste("GPP (g C m"^-2, " d"^-1, ")" ) ),
       x = "Day of year") +
  facet_wrap(~sitename)

ggsave("./fig/meanseasonalcycle_by_site_cal_fromKoen.pdf", width = 15, height = 30)


