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
    
    cwd = os.getcwd()
    eps0 = model.model_param_values
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    fitit = 1
    
    
    
    ##calculate the scaling based upon what would actually happen: if it's already at 0.01, and goes more negative, now it effectively shows a deps there of 0, and not some arbitrarily large number
    neweps_effective = eps0 + deps
    #deps_actual = neweps_effective - eps0
    neweps_effective[neweps_effective < 0.01] = 0.01
    deps_effective = neweps_effective - eps0
    
    
    factor = np.linalg.norm(deps_effective)/np.linalg.norm(eps0)
    max_step = np.max(np.abs(deps_effective))
    
        
    if factor > max_step_factor:
        deps = (deps*max_step_factor) / max_step
        print "Scaling down to %f by maximum step" % max_step_factor
        fitit = max_step_factor/max_step
    
    neweps = eps0 + deps
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
    
    eplot.plot_epsilons_bin(deps,"d-epsilon",model)
    eplot.plot_epsilons(deps,"d-epsilon",model)
    
    


        
        
        
        
        
        
        
        
        
        
        
    

