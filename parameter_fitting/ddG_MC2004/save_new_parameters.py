""" Save the new parameters from the Newton solution


"""

import numpy as np
import os



def save(model,soln_index):
    """ """

    eps0 = model.epsilons
    deps = np.loadtxt("xp_%d.dat" % soln_index)


    
