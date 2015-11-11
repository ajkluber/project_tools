''' Save the new parameters from the Newton solution


'''

import numpy as np
import os

def save(model,fitopts,soln_index,nonnative=False):
    ''' Save new parameters '''

    # Only the fitting_params (a subset of model_params) are
    # being updated.

    # Use new method of model_parameters
    eps_p_0 = model.model_param_values[model.fitting_params]
    deps_p = np.loadtxt("xp_%d.dat" % soln_index)

    # Scale model parameters by a constant such that the ratio 
    # of the norm of the perturbation to the norm of the parameters
    # is 0.2.
    if os.path.exists("desired_ratio"):
        desired_ratio = np.loadtxt("desired_ratio")
    else:
        desired_ratio = 1.
    print " scaling step size to %.3f" % desired_ratio

    nrm_model_params = np.linalg.norm(eps_p_0)
    ratio = np.linalg.norm(deps_p)/nrm_model_params
    if ratio > desired_ratio:
        alpha = desired_ratio/ratio
    else:
        alpha = 1.
    open("scaling_alpha","w").write("%.4f" % alpha)
    neweps_p = eps_p_0 + alpha*deps_p       
    
    # For the non-native interactions set the neweps_p to be
    # equal to the delta.
    #for i in range(model.n_contacts):
    #    #
    #    pass 

    # Update parameters
    if nonnative:
        update_model_param_values_nonnative(model, neweps_p)
    else:
        model.update_model_param_values(neweps_p)

    cwd = os.getcwd()
    relpath = cwd.split("%s/" % model.path)[1]
    open("pairwise_params","w").write(model.pairwise_param_file_string)
    open("model_params","w").write(model.model_param_file_string)
    model.pairwise_params_file_location = "%s/pairwise_params" % relpath
    model.model_params_file_location = "%s/model_params" % relpath



def update_model_param_values_nonnative(model,new_model_param_values):
    """ If parameter changed sign, change the pairwise interaction type """
    # Switching between different interaction function types
    potential_type_switch = {4:9,9:4}

    # Loop over fitting_params only 
    for i in range(model.n_fitting_params):
        p_idx = model.fitting_params[i]
        p_pairs = model.model_param_interactions[p_idx]
        for n in range(len(p_pairs)):
            if new_model_param_values[i] < 0.:
                # If parameter changes sign then the pairwise_type is flipped.
                if model.pairwise_type[p_pairs[n]] == 1:
                    model.pairwise_type[p_pairs[n]] = 2
                else:
                    model.pairwise_type[p_pairs[n]] = potential_type_switch[model.pairwise_type[p_pairs[n]]]

            # Model parameters are always positive
            model.model_param_values[p_idx] = abs(new_model_param_values[i])   

    # Refresh everything that depends on model parameters
    model._set_nonbonded_interactions()
