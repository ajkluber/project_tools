""" Analysis scripts for analyzing coarse-grain simulations

Description:

    The analysis submodule holds a collection of analysis code, some which can
be used as stand-alone utilities.


Submodules:
    
crunch_coordinates -- Functions that submit PBS jobs to calculate everything.

wham -- Preps and submits PBS jobs to perform WHAM (1) calculations.

whamdata -- Holds a class specialized to data input file data for Cecilia's wham code.

References:
(1) Kumar, S.; Bouzida, D.; Swendsen, R. H.; Kollman, P. A.; Rosenberg, J. The 
Weighted Histogram Analysis Method for Free-energy Calculations on Biomolecules
I. The Method. J. Comp. Chem. 1992, 13, 1011-1021
"""

import crunch_coordinates
import constant_temp
import wham
import smog_AA_constant_temp
import smog_AA_crunch_coordinates
