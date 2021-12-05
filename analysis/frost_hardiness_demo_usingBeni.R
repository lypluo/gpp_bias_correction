#----------------------------------------------------
# run a small example of the frost hardiness function-->using the Beni's hardening function
#----------------------------------------------------
library(tidyverse)
library(GenSA)
source("R/model_hardening_byBeni.R")

# read in data
df <- readRDS("data/model_data.rds") %>%
  mutate(
    year = format(date, "%Y")
  ) %>%
  filter(
    sitename == "US-Ha1",
    date > "2012-01-01"
  )

#set initial parameters for the model
par <- c("a" = 0, "b" = 0.5, "c" = 50, "d" = 0.1, "e" = 1)

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
      scaling_factor = scaling_factor
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
      gpp_mod * scaling_factor
    ),
    colour = "blue"
  ) +
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
