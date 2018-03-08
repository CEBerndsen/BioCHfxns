#Function to predict Michaelis-Menten Kinetics
library(tidyverse)
library(broom)

  MM_simplefit_plot <- function(df){
  #df should be a two column data frame with columns names 'rate' and 'conc'
  df <- df
  #Establish some values for fitting
  Kstart <- max(df$conc)/2
  Vstart <- max(df$rate)
  
  #Fit the data
  fit <- nls(rate ~ (V * conc)/(K + conc), data = df, start = list(K = Kstart, V = Vstart))
  
  #Tidy the table
  fitstats <- tidy(fit)
  
  #make the fit data
  mmfitdata <- data.frame(conc = seq(0, max(df$conc)*1.2, length.out = 200))
  mmfitdata$rate <- predict(fit, newdata = mmfitdata)
  
  ggplot() +
      geom_point(data = df, aes(x = conc, y = rate), color = "red") + 
      geom_line(data = mmfitdata, aes(x = conc, y = rate), colour = "black")
    
}


df <- data_frame(conc = seq(0, 100, by = 10),
                 rate = c(0, 3, 6, 9, 12, 14, 15, 16, 17, 17.5, 18))


