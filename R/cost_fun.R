##cost functions for different model structures

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

#remove the parameter e
cost_new4 <- function(data,par_new4) {
  scaling_factor <- data %>%
    # group_by(sitename) %>%
    do({scaling_factor <- model_hardening_new4(.,par_new4)
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
