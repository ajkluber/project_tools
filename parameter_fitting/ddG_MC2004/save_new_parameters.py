""" Save the new parameters from the Newton solution


"""

import numpy as np
import os


def save(model,soln_index):
    """ Save new parameters """
    ## See model_builder/models/pairwise_potentials for codes
    potential_type_switch = {1:2,2:3,3:2}

    ## Use new method of model_parameters
    eps_p_0 = model.model_param_values
    deps_p = np.loadtxt("xp_%d.dat" % soln_index)

    ## Scale model parameters by a constant such that the ratio 
    ## of the norm of the perturbation to the norm of the parameters
    ## is 0.1.
    nrm_model_params = np.linalg.norm(eps_p_0)
    ratio = np.linalg.norm(deps_p)/nrm_model_params
    alpha = 0.1/ratio
    neweps_p = eps_p_0 + alpha*deps_p       
    
    ## NEED TO DO:
    ##  - Determine which contacts switched from attractive to repulsive.
    ##    and change the function type in model accordingly

    ## Update parameters
    model.model_param_values = neweps
    model.generate_topology()   
    model.contact_params = "%s/NewBeadBead.dat" % cwd

    ## vvvv DEPRECATED vvvv
    #if model.fitting_allowswitch ==False:
    #    neweps[neweps < 0.01] = 0.01
    #else:
    #    index = np.where(neweps<0.)
    #    for i in range(len(index)):
    #        # If epsilon becomes negative, change interaction type. Keep epsilon positive
    #        model.LJtype[index[i]] = -model.LJtype[index[i]] 
    #        neweps[index[i]] = abs(neweps[index[i]])
    #model.contact_epsilons = neweps
    #model.generate_topology()
    #open("NewBeadBead.dat","w").write(model.beadbead)
    #model.contact_params = "%s/NewBeadBead.dat" % cwd
