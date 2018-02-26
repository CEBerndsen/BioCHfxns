# BioCHfxns
Functions for biochemistry

A series of R functions for use in biochemical and biophysical experiments. 

# List of functions
CD_predict: A tool that uses the work of Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629 to predict circular dichroism spectra. Inputs are fraction alpha helix (a), fraction beta strand/sheet (b) and fraction random coil (c) which must sum to 1 (a + b + c == 1). if guides are set to TRUE then the canonical spectra for an alha helical protein, a beta sheet protein or a random coil protein are overlaid on the predicted spectrum.
