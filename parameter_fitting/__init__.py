""" Parameter fitting using Newton's method

Description:

    This submodule performs parameter fitting to match simulation features to
experimental (or designed) data. This parameter fitting is an extension of the
Matysiak, Clementi 2004 (see reference (1)). 


Submodules:

  newton 
      Solves for new parameters by Newton's method: by inverting the 
    Jacobian.

  ddG_MC2004
      Calculates Jacobian and feature vector for Delta Delta G's to match
    experimental data from phi-value analysis.


References:

1. Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.

2. J.E. Dennis; Robert B. Schnabel. "Numerical Methods for Unconstrained
Optimization and Nonlinear Equations". SIAM. 1996.

3. Jorge Nocedal; Stephen Wright. "Numerical Optimization". Springer. 2000


"""

import newton
import ddG_MC2004
import FRET
import RMSF
