""" Save the new parameters from the Newton solution


"""

import numpy as np
import os


def save(model,soln_index):
    """ Save new parameters """
    
    cwd = os.getcwd()
    eps0 = model.pairwise_strengths
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    diff_factor = np.absolute(deps / eps0)
    
    model.pairwise_strengths = neweps
    model.generate_topology()
    open("NewBeadBead.dat","w").write(model.beadbead)

    model.contact_params = "%s/NewBeadBead.dat" % cwd
