#Function to predict inhibition type
library(tidyverse)
library(broom)

#predict function

EnzInh_predict <- function(df, Vs = max(df$rate)*1.5, Ks = max(df$conc)/2, details = FALSE) {
  #df should be a three column data frame with column names 'rate' and 'conc' and at least two 'Inhibitor concentrations (I)'
  df <- df

  #Fit the data
  fitcomp <- nls(rate ~ (V * conc)/(K*(1 + (I/Ki)) + conc), data = df, start = list(K = Ks, V = Vs, Ki = Ks))
  fitnoncomp <- nls(rate ~ (V/(1 + (I/Ki)) * conc)/(K + conc), data = df, start = list(K = Ks, V = Vs, Ki = Ks))
  fituncomp <- nls(rate ~ (V/(1 + (I/Ki)) * conc)/(K/(1 + (I/Ki)) + conc), data = df, start = list(K = Ks, V = Vs, Ki = Ks))
  
  if(details == TRUE){
    dat <- list(fitcomp, fitnoncomp, fituncomp)
    names(dat) = c("competitive", "noncompetitive", "uncompetitive")
    return(dat)
     }
  else{
    #tidy tables and annotate mechanism
    fitcomp <- fitcomp %>% tidy(.) %>% mutate(mechanism = "competitive") 
    fitnoncomp <- fitnoncomp %>% tidy(.) %>% mutate(mechanism = "noncompetitive")
    fituncomp <- fituncomp %>% tidy(.) %>% mutate(mechanism = "uncompetitive")
    #make comparison table
    mechtbl <- bind_rows(fitcomp, fitnoncomp, fituncomp) %>% select(-statistic, -p.value)
    return(mechtbl)
  }
}

