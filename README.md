# BioCHfxns
Functions for biochemistry

A series of R functions for use in biochemical and biophysical experiments. 

# List of functions
CD_predict: A tool that uses the work of Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629 to predict circular dichroism spectra. Inputs are fraction alpha helix (a), fraction beta strand/sheet (b) and fraction random coil (c) which must sum to 1 (a + b + c == 1). if guides are set to TRUE then the canonical spectra for an alha helical protein, a beta sheet protein or a random coil protein are overlaid on the predicted spectrum.

MM_simplefit: A function that takes a data frame with columns conc and rate. Returns a table showing the Vmax and Km values plus error

MM_simplefit_plot: A function that takes a data frame with columns conc and rate. Returns a ggplot2 plot showing the data plus the fit as a line. Plot object can be modifed by assigning the MM_simplefit_plot object to a variable and then using the '+' to add on theme or other geoms.
