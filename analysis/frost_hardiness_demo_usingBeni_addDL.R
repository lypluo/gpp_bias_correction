#----------------------------------------------------
# run a small example of the frost hardiness function-->using the Beni's hardening function
#----------------------------------------------------
library(tidyverse)
library(GenSA)
#compare two source of stress functions:
source("R/model_hardening_byBeni.R")
source("R/model_hardening_byBeni_addDL.R")

# read in data
df <- readRDS("data/model_data.rds") %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  filter(
    sitename == "US-Ha1",
    date > "2012-01-01"
  )

#add day length into the datasets:
library(geosphere) #calculate the day length
library(sirad)  #calculate the day of the year
library(lubridate)
#load the lat and long of the example site:
df_sites_modis_era <-read_rds(paste0("./data/df_sites_modis_era.csv"))
siteinfo<-df_sites_modis_era %>% filter(sitename=="US-Ha1")
#
df$doy<-dayOfYear(df$date)
df$DL<-daylength(siteinfo$lat,df$doy)

#set initial parameters for the model
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1,"f"=6)

#function without adding DL
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening(
      .,
      par
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_noDL = scaling_factor
    )
  })

df <- left_join(df, scaling_factors)

#function with adding DL
scaling_factors <- df %>%
  group_by(sitename, year) %>%
  do({
    scaling_factor <- model_hardening_new(
      .,
      par
    )

    data.frame(
      sitename = .$sitename,
      date = .$date,
      scaling_factor_DL = scaling_factor
    )
  })
df <- left_join(df, scaling_factors)

# exploratory plot
p <- ggplot(df) +
  geom_line(
    aes(
      date,
      gpp_mod
    ),
    colour = "red"
  ) +
  geom_line(
    aes(
      date,
      gpp_mod * scaling_factor_noDL
    ),
    colour = "blue"
  ) +
  geom_line(
    aes(
      date,
      gpp_mod * scaling_factor_DL
    ),
    colour = "orange"
  )+
  geom_point(
    aes(
      date,
      gpp
    )
  ) +
  labs(
    x = "",
    y = "GPP",
    title = "US-NR1"
  ) +
  theme_minimal()

# ggsave("niwot.png", width = 7, height = 4)
