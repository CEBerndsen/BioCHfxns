#Function to predice Michaelis-Menten Kinetics
library(tidyverse)
library(broom)

MM_simplefit <- function(df, guides = TRUE, plot = TRUE){
  #df should be a two column data frame with columns names 'rate' and 'conc'
  
  #Establish some values for fitting
  Kstart <- max(df$conc)/2
  Vstart <- max(df$rate)
  
  #Fit the data
  fit <- nls(rate ~ (V * conc)/(K + conc), data = df, start = list(K = Kstart, V = Vstart))
  
  #Tidy the table
  fitstats <<- tidy(fit)
  
  #make the fit data
  mmfitdata <- data.frame(conc = seq(0, max(df$conc)*1.2, length.out = 200))
  mmfitdata$Rate <- predict(fit, newdata = mmfitdata)
  
  if(plot == TRUE){
  #make a plot
    if(guides == TRUE) {
    ggplot() +
      geom_point(data = df, aes(x = conc, y = rate)) + 
      geom_line(data = mmfitdata, aes(x = Conc, y = Rate*1000), colour = "black") +
      geom_segment(aes(x = 0, y = (0.7021/2), xend = 210, yend = (0.7021/2)), color = "red") +
      geom_segment(aes(x = 210, y = 0, xend = 210, yend = (0.7021/2)), color = "red") +
      geom_hline(yintercept = 0.7021, color = "purple") +
      annotate("text", x = 2500, y = 0.8, label = "~V[max]", color = "purple", parse = TRUE, size = 8) +
      annotate("text", x = 600, y = 0.1, label = "~K[M]", color = "red", parse = TRUE, size = 8) +
      geom_segment(aes(x = 450, y = 0.07, xend = 250, yend = 0.01), color = "red", arrow = arrow(length = unit(0.1, "inches")))
  } else
  {
      ggplot() +
        geom_point(data = df, aes(x = conc, y = rate), color = "red") + 
        geom_line(data = mmfitdata, aes(x = Conc, y = Rate*1000), colour = "black")
    }
  } else
  {
    return(fitstats)
    }
}
