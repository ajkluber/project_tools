""" Save the new parameters from the Newton solution


"""

import numpy as np
import os


def save(model,soln_index):
    """ Save new parameters """
    
    cwd = os.getcwd()
    eps0 = model.pairwise_strengths
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    neweps = eps0 + deps
    neweps[neweps < 0.01] = 0.01

    ## Take choice and rescale it so that the smallest new 
    ## parameter is epsmin.
    #Alphas = -(eps0 - epsmin)/deps
    #alpha = min(Alphas[Alphas > 0])
    #if alpha > 1.0:
    #    neweps = eps0 + deps
    #else:
    #    neweps = eps0 + alpha*deps

    ## Decide on whether to keep new repulsive contacts
    #if model.fitting_allowswitch == "True":
    #   pass
    #   ## Make negative epsilons positive and change the
    #   ## the interaction to repulsive.
    #elif model.fitting_allowswitch == "False":
    #   pass
    #   neweps[neweps < 0.01] = 0.01
    #else:
    #   print "ERROR!"
    #if model.fitting_allowswitch ==False:
    #    neweps[neweps < 0.01] = 0.01
    #else:
    #    index = np.where(neweps<0.)
    #    for i in range(len(index)):
    #        # If epsilon becomes negative, revert contact type. Keep epsilon positive
    #        model.LJtype[index[i]] = -model.LJtype[index[i]] 
    #        neweps[index[i]] = abs(neweps[index[i]])
        
    model.pairwise_strengths = neweps
    model.generate_topology()
    open("NewBeadBead.dat","w").write(model.beadbead)

    model.contact_params = "%s/NewBeadBead.dat" % cwd
