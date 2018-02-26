#Load packages
library(tidyverse)


predict_CD <- function(a, b, c, guides = TRUE) {
  #check to make sure fractions add to 1          
  stopifnot(a <= 1, b <= 1, c <= 1, a + b + c == 1)
  
  #Generate the wavelength values
  CDdat <- data.frame(lambda = seq(190, 250, by = 0.2))
  
  #Generate the basis set from Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629
  CDdat <- CDdat %>% mutate(helix = 1*10^8 * (2230060.04151075*lambda^0 +
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
  
  #Predict spectrum based on user input  
  CDdat <- CDdat %>% mutate(prediction = a*helix + b*beta + c*coil)
  
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



# finding missing parenthesis ---------------------------------------------

if_else(input$protein == 'lyso',
        CDdat %>% 
          mutate(prediction = 0.7*helix + 0.07*beta + 0.17*coil) %>%
          ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
          geom_jitter(alpha = 0.4) +
          scale_x_continuous(breaks = seq(190, 250, by = 5)) +
          labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for lysozyme", caption = "70% helix, 7% Beta strand, 17% random coil") +
          geom_hline(yintercept = 0) +
          ylim(-50, 75) +
          theme_classic() +
          theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")),
        if_else(input$protein == 'ub',
                CDdat %>% 
                  mutate(prediction = 0.17*helix + 0.40*beta + 0.43*coil) %>%
                  ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
                  geom_jitter(alpha = 0.4) +
                  scale_x_continuous(breaks = seq(190, 250, by = 5)) +
                  labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for ubiquitin", caption = "17% helix, 40% Beta strand, 43% random coil") +
                  geom_hline(yintercept = 0) +
                  ylim(-50, 75) +
                  theme_classic() +
                  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")),
                if_else(input$protein == 'bst',
                        CDdat %>% 
                          mutate(prediction = 0.92*helix + 0*beta + 0.08*coil) %>%
                          ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
                          geom_jitter(alpha = 0.4) +
                          scale_x_continuous(breaks = seq(190, 250, by = 5)) +
                          labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for BST2", caption = "92% helix, 0% Beta strand, 8% random coil") +
                          geom_hline(yintercept = 0) +
                          ylim(-50, 75) +
                          theme_classic() +
                          theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")),
                        if_else(input$protein == 'hemo',
                                CDdat %>% 
                                  mutate(prediction = 0.74*helix + 0*beta + 0.17*coil) %>%
                                  ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
                                  geom_jitter(alpha = 0.4) +
                                  scale_x_continuous(breaks = seq(190, 250, by = 5)) +
                                  labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for hemoglobin", caption = "74% helix, 0% Beta strand, 17% random coil") +
                                  geom_hline(yintercept = 0) +
                                  ylim(-50, 75) +
                                  theme_classic() +
                                  theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")),
                                if_else(input$protein == 'ab',
                                        CDdat %>% 
                                          mutate(prediction = 0.05*helix + 0.46*beta + 0.36*coil) %>%
                                          ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
                                          geom_jitter(alpha = 0.4) +
                                          scale_x_continuous(breaks = seq(190, 250, by = 5)) +
                                          labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for an antibody", caption = "5% helix, 46% Beta strand, 36% random coil") +
                                          geom_hline(yintercept = 0) +
                                          ylim(-50, 75) +
                                          theme_classic() +
                                          theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")),
                                        ggplot() +
                                          geom_line(data = CDdat, aes(x = lambda, y = helix/1000), color = "red") +
                                          geom_line(data = CDdat, aes(x = lambda, y = beta/1000), color = "green") +
                                          geom_line(data = CDdat, aes(x = lambda, y = coil/1000), color = "purple") +
                                          geom_hline(yintercept = 0) +
                                          scale_x_continuous(breaks = seq(190, 250, by = 5)) +
                                          labs(x = "wavelength (nm)", y = "Ellipticity", title = "CD spectrum") +
                                          annotate("text", x = 235, y = 25, label = "Helix", color = "red", size = 5) +
                                          annotate("text", x = 235, y = 20, label = "Beta Sheet", color = "green", size = 5) +
                                          annotate("text", x = 235, y = 15, label = "Random Coil", color = "purple", size = 5) +
                                          ylim(-50, 75) +
                                          theme_classic() +
                                          theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold")))))))


# make the data set -------------------------------------------------------

#Generate the wavelength values
CDdat <- data.frame(lambda = seq(190, 250, by = 0.2))

#Generate the basis set from Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629
CDdat <- CDdat %>% mutate(helix = 1*10^8 * (2230060.04151075*lambda^0 +
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

#Predict spectrum based on user input  
CDdat <- CDdat %>% 
  mutate(lyso = 0.7*helix + 0.07*beta + 0.17*coil) %>%
  mutate(ub = 0.17*helix + 0.40*beta + 0.43*coil) %>%
  mutate(bst = 0.92*helix + 0*beta + 0.08*coil) %>%
  mutate(hemo = 0.74*helix + 0*beta + 0.17*coil) %>%
  mutate(ab = 0.05*helix + 0.46*beta + 0.36*coil)



# naming code -------------------------------------------------------------

x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for an antibody", caption = "5% helix, 46% Beta strand, 36% random coil"

labels <- ifelse(input$protein == 'lyso', x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for lysozyme", caption = "70% helix, 7% Beta strand, 17% random coil",
                  ifelse(input$protein == 'ub', x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for ubiquitin", caption = "17% helix, 40% Beta strand, 43% random coil",
                          ifelse(input$protein == 'bst', x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for BST2", caption = "92% helix, 0% Beta strand, 8% random coil",
                                  ifelse(input$protein == 'ab', x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for an antibody", caption = "5% helix, 46% Beta strand, 36% random coil",
                                         ifelse(input$protein == 'hemo', x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum for hemoglobin", caption = "74% helix, 0% Beta strand, 17% random coil", x = "wavelength (nm)", y = "Ellipticity")))))
  


label <- ifelse(input$protein == 'lyso', "70% helix, 7% Beta strand, 17% random coil",
                ifelse(input$protein == 'ub', "17% helix, 40% Beta strand, 43% random coil",
                       ifelse(input$protein == 'bst', "92% helix, 0% Beta strand, 8% random coil",
                              ifelse(input$protein == 'ab', "5% helix, 46% Beta strand, 36% random coil",
                                     ifelse(input$protein == 'hemo', "74% helix, 0% Beta strand, 17% random coil", "N/A")))))

output$text <- renderText({
  paste("Input text is:", label)