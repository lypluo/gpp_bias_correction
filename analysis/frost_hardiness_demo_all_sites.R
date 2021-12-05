# run a small example of the frost hardiness
# function

library(tidyverse)
library(GenSA)
source("R/frost_hardiness.R")

# read in data
df <- readRDS("data/model_data.rds") %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  na.omit()

# note: series need to be complete, NAs are not
# allowed, but this is a quick fix

# first order model for forst hardiness
# after Leinone et al. 1995 and simplified by
# Hanninen and Kramer 2007
# literature values of the
# frost hardiness function

par <- c(
  a = 1.5,
  b = -21.5,
  T8 = 11.3,
  t = 5,
  base_t = -25
)

scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- frost_hardiness(
      .$temp,
      par
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor = scaling_factor
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
    geom_line(aes(date,gpp_mod * scaling_factor,col="cal-mod")) +
    labs(x = "",y = "GPP") +
    # facet_grid(sitename ~ .) +
    scale_color_manual(values = c("ori-mod"="red","cal-mod"="steelblue2","obs"="grey"))+
    annotate(geom = "text",x=label_x,y=label_y,label=site_run)+
    theme_minimal()
  p_all[[i]]<-p
}
names(p_all)<-site_names

# ggsave("~/Desktop/time_series.pdf", height = 48, width = 8)
