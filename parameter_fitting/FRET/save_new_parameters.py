""" Save the new parameters from the Newton solution


"""

import numpy as np
import os
import matplotlib.pyplot as plt

def save(model,fitopts,soln_index):
    """ Save new parameters """
    import analysis_scripts.plot_depsilon_native as eplot
    
    if "max_step_factor" in fitopts:
        max_step_factor = float(fitopts["max_step_factor"])
    else:
        max_step_factor = 0.3
    
    if "prevent_zero" in fitopts:
        prevent_zero = fitopts["prevent_zero"]
    else:
        prevent_zero = True
    
    if model.using_sbm_gmx == None or model.using_sbm_gmx == False:
        gauss = False
    elif model.using_sbm_gmx == True:
        gauss = True
    else:
        gauss = False
    
    cwd = os.getcwd()
    eps0 = model.model_param_values
    temp_deps = np.loadtxt("xp_%d.dat" % soln_index)
    fitit = 1
    ##read in the deps and assign them to the proper matrix for telegraphing to eps0 correctly
    deps = np.zeros(model.n_model_param)
    for i in range(model.n_fitting_params):
        param_idx = model.fitting_params[i]
        deps[param_idx] = temp_deps[i]
    
    
    ##calculate the scaling based upon what would actually happen: if it's already at 0.01, and goes more negative, now it effectively shows a deps there of 0, and not some arbitrarily large number
    neweps_effective = eps0 + deps
    #deps_actual = neweps_effective - eps0
    if prevent_zero:
        neweps_effective[neweps_effective < 0.01] = 0.01
    deps_effective = neweps_effective - eps0
    
    
    factor = np.linalg.norm(deps_effective)/np.linalg.norm(eps0)
    max_step = np.max(np.abs(deps_effective))
    
        
    if factor > max_step_factor:
        deps = (deps*max_step_factor) / max_step
        print "Scaling down to %f by maximum step" % max_step_factor
        fitit = max_step_factor/max_step
    
    neweps = eps0 + deps
    if prevent_zero:
        neweps[neweps < 0.01] = 0.01
    
    model.update_model_param_values(neweps)
    
    plt.figure()
    plt.hist(neweps, 50, alpha=0.75)
    plt.xlabel("epsilon",fontsize=20)
    plt.ylabel("number",fontsize=20)
    plt.title("spread of epsilons", fontsize=20)
    plt.savefig("eps_spread.png")
    
    open("%s/pairwise_params" % cwd,"w").write(model.pairwise_param_file_string)
    open("%s/model_params" % cwd,"w").write(model.model_param_file_string)
    open("%s/fitting_scale" % cwd,"w").write("%f\n%f"%(fitit, np.sum(neweps)/np.shape(neweps)[0]))
    model.pairwise_params_file_location = "%s/pairwise_params" % cwd
    model.model_params_file_location = "%s/model_params" % cwd
    
    if gauss:
        eplot.plot_epsilons_bin(deps[np.arange(1,len(deps),2)],"d-epsilon",model,gauss=gauss)
        eplot.plot_epsilons(deps[np.arange(1,len(deps),2)],"d-epsilon",model,gauss=gauss)
    else:
        eplot.plot_epsilons_bin(deps,"d-epsilon",model)
        eplot.plot_epsilons(deps,"d-epsilon",model)
    
    


        
        
        
        
        
        
        
        
        
        
        
    

