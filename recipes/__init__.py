""" Submodule containing recipes for running specific algorithms.

Description:

    This submodule is called 'recipes' because it holds scripts to
automatically execute a set of steps to implement a specific algorithm. The
purpose of this submodule is ensure that procedures are easily reproduceable
and extendable to additional systems.


Available Recipes:

MatysiakClementi2004
    Implements procedure from reference (1) to get a heterogeneous Go-model.
    
FRETfit (coming soon)
    Implements algorithm to match a heterogeneous Go-model to FRET data.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""
import recipe_manager
import MatysiakClementi2004
import Clementi2000
import FRETFit 
