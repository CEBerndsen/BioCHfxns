#Function to predict Michaelis-Menten Kinetics
library(tidyverse)
library(broom)


MM_simplefit <- function(df){
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
  
  return(fitstats)
  
}