""" Save the new parameters from the Newton solution


"""

import numpy as np
import os


def save(model,soln_index):
    """ Save new parameters """
    
    cwd = os.getcwd()
    eps0 = model.contact_epsilons
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    epsmin = 0.01*np.ones(len(eps0),float)

    ## Take choice and rescale it so that the smallest new 
    ## parameter is epsmin.
    Alphas = -(eps0 - epsmin)/deps
    alpha = min(Alphas[Alphas > 0])

    neweps = eps0 + alpha*deps

    model.contact_epsilons = neweps
    model.generate_topology()
    open("NewBeadBead.dat","w").write(model.beadbead)

    model.contact_params = "%s/NewBeadBead.dat" % cwd
