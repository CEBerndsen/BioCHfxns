#Load packages
library(tidyverse)
library(readr)
library(broom)

predict_CD <- function(df, guides = TRUE) {
  #Load the dataset, data should be in the form of wavelength in one column and ellipticity in the second column
  CDdat <- read.csv(df)
  colnames(CDdat) <- c("lambda", "data")
  
  #Generate the data frame 
  CDdf <- CDdat %>% 
    #remove data outside good prediction limits
    filter(lambda <= 250) %>%
    filter(lambda > 190) %>%
    #Generate the basis set basis set from Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629
    mutate(helix = 1*10^8 * (2230060.04151075*lambda^0 +
                                                -100548.516559741*lambda^1 +
                                                2037.18080475746*lambda^2 + 
                                                -24.4244919907991*lambda^3 +
                                                0.19190243015954*lambda^4 +
                                                -0.00103245782924168*lambda^5 +
                                                0.00000385211889091252*lambda^6 +
                                                -9.84175959744622E-09*lambda^7 +
                                                1.64786777298595E-11*lambda^8 +
                                                -1.63282751503442E-14*lambda^9 +
                                                7.27089674019501E-18*lambda^10)) %>% 
    mutate(beta = 1*10^8 * (-677807.330017282*lambda^0 +
                              30975.2707887604*lambda^1 +
                              -636.143263740698*lambda^2 + 
                              7.73164864362657*lambda^3 +
                              -6.15861633716145E-02*lambda^4 +
                              3.35943314432255E-04*lambda^5 +
                              -1.27092416556044E-06*lambda^6 +
                              3.29272089581372E-09*lambda^7 +
                              -5.59118151750062E-12*lambda^8 +
                              5.61899389095424E-15*lambda^9 +
                              -2.53794637234403E-18*lambda^10)) %>%
    mutate(coil = 1*10^8 * (-580939.072386969*lambda^0 +
                              25845.2673351998*lambda^1 +
                              -516.713088253122*lambda^2 + 
                              6.1134023680003*lambda^3 +
                              -4.74021175198809E-02*lambda^4 +
                              2.51692531821056E-04*lambda^5 +
                              -9.26824208397782E-07*lambda^6 +
                              2.33714935193268E-09*lambda^7 +
                              -3.86247107852678E-12*lambda^8 +
                              3.77764956561175E-15*lambda^9 +
                              -1.6603998403172E-18*lambda^10))
  
  #fit the data spectrum based on user input  
  CDmod <- nls(data ~ a * helix + b * beta + c * coil, data = CDdf, start = list(a = 0.2, b = 0.5, c = 0.1))
  
  CDdf <- CDdf %>% mutate(prediction = 3.29e-7*helix + 3.875e-8*beta + 2.157e-7*coil)
  
  
  #Plot the prediction and save as an object
  if(guides == FALSE) {
    ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
      geom_jitter(alpha = 0.4) +
      scale_x_continuous(breaks = seq(190, 250, by = 5)) +
      labs(x = "wavelength (nm)", y = "Ellipticity") +
      geom_hline(yintercept = 0) +
      theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold"))
  }
  else {
    ggplot() +
      geom_jitter(data = CDdat, aes(x = lambda, y = prediction/1000), fill = "red", alpha = 0.4) +
      geom_line(data = CDdat, aes(x = lambda, y = helix/1000), color = "red") +
      geom_line(data = CDdat, aes(x = lambda, y = beta/1000), color = "green") +
      geom_line(data = CDdat, aes(x = lambda, y = coil/1000), color = "purple") +
      scale_x_continuous(breaks = seq(190, 250, by = 5)) +
      labs(x = "wavelength (nm)", y = "Ellipticity") +
      theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold"))
  }
  
  
  
  
}

