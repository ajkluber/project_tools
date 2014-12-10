''' Save the new parameters from the Newton solution


'''

import numpy as np
import os


def save(model,soln_index):
    ''' Save new parameters '''
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
    
    ## For the non-native interactions set the neweps_p to be
    ## equal to the delta.
    for i in range(model.n_contacts):

    ## Update parameters
    model.update_model_param_values(neweps_p)

    open("%s/pairwise_params" % cwd,"w").write(self.pairwise_param_file)
    open("%s/model_params" % cwd,"w").write(self.model_param_file)
    model.pairwise_param_file = "%s/pairwise_params" % cwd
    model.model_param_file = "%s/model_params" % cwd

