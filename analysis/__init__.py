""" Analysis scripts for analyzing coarse-grain simulations

Description:

    The analysis submodule holds a collection of analysis code, some which can
be used as stand-alone utilities.


Submodules:
    
contacts -- Code for computing residue-residue contacts. 

crunch_coordinates -- Functions that submit PBS jobs to calculate everything.

plot -- Submodule to make various plots, like contact maps and free energies.    

Tf_loop -- Routines that are used in recipes. Not really stand-alone.

wham -- Preps and submits PBS jobs to perform WHAM (1) calculations.

whamdata -- Holds a class specialized to data input file data for Cecilia's wham code.


Other excutables:

WHAM
    Binary executable of Cecilia's fortran90 code to perform Weighted Histogram
Analysis Method (WHAM) over a set of temperatures. See reference (1) for details
about the WHAM method.


References:
(1) Kumar, S.; Bouzida, D.; Swendsen, R. H.; Kollman, P. A.; Rosenberg, J. The 
Weighted Histogram Analysis Method for Free-energy Calculations on Biomolecules
I. The Method. J. Comp. Chem. 1992, 13, 1011-1021
"""

import contacts 
import crunch_coordinates
import plot
import constant_temp
import wham
import TPT
