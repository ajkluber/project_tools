''' Save the new parameters from the Newton solution


'''

import numpy as np
import os
import shutil

def save(model,fitopts,soln_index,nonnative=False):
    ''' Save new parameters '''

    # Only the fitting_params (a subset of model_params) are
    # being updated.

    # Use new method of model_parameters
    eps_p_0 = model.long_model_param_values[model.long_fitting_params]
    deps_p = np.loadtxt("xp_%d.dat" % soln_index)

    # Scale model parameters by a constant such that the ratio 
    # of the norm of the perturbation to the norm of the parameters
    # is equal to "ratio"

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
    open("smog_pairs_l.top","w").write(model.long_pairs_file_string)
    open("smog_pairs_long","w").write(model.long_pairs_file_string)
    open("smog_bonds_rep.top","w").write(model.long_bonds_rep_string)
    shutil.move('../../smog_files/smog_pairs_long','../../smog_files/smog_pairs_long_{0}'.format(int(fitopts['iteration'])))
    shutil.move('../../smog_files/smog_pairs_l.top','../../smog_files/smog_pairs_l_top_{0}'.format(int(fitopts['iteration'])))
    shutil.move('../../smog_files/smog_bonds_rep.top','../../smog_files/smog_bonds_rep_top_{0}'.format(int(fitopts['iteration'])))
    shutil.copy('smog_pairs_l.top','../../smog_files/smog_pairs_l.top')
    shutil.copy('smog_pairs_long','../../smog_files/smog_pairs_long')
    shutil.copy('smog_bonds_rep.top','../../smog_files/smog_bonds_rep.top')
#    model.long_pairs_file_location = "%s/smog_pairs_long" % relpath

def update_model_param_values_nonnative(model,new_model_param_values):
    """ If parameter changed sign, change the pairwise interaction type """
    # Switching between different interaction function types
    potential_type_switch = {5:9,9:5}

    # Loop over fitting_params only 
    for i in range(model.n_long_fitting_params):
        p_idx = model.long_fitting_params[i]
        if new_model_param_values[i] < 0.:
                # If parameter changes sign then the pairwise_type is flipped.
            model.long_pairwise_type[p_idx] = potential_type_switch[model.long_pairwise_type[p_idx]]
            # Model parameters are always positive
        else:
            pass

        model.long_model_param_values[p_idx] = abs(new_model_param_values[i])

    # Refresh everything that depends on model parameters
    model._determine_tabled_interactions()
    model._set_nonbonded_interactions()
