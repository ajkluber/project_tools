""" Save the new parameters from the Newton solution


"""

import numpy as np
import os


def save(model,fitopts,soln_index):
    """ Save new parameters """
    
    cwd = os.getcwd()
    eps0 = model.contact_epsilons
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    epsmin = 0.01*np.ones(len(eps0),float)
    neweps = eps0 + deps

    ## Take choice and rescale it so that the smallest new 
    ## parameter is epsmin.
    #Alphas = -(eps0 - epsmin)/deps
    #alpha = min(Alphas[Alphas > 0])
    #if alpha > 1.0:
    #    neweps = eps0 + deps
    #else:
    #    neweps = eps0 + alpha*deps

    neweps[(neweps < 0.01)] = 0.01

    model.contact_epsilons = neweps
    model.generate_topology()
    open("NewBeadBead.dat","w").write(model.beadbead)

    model.contact_params = "%s/NewBeadBead.dat" % cwd
