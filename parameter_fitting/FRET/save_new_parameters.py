""" Save the new parameters from the Newton solution


"""

import numpy as np
import os
import matplotlib.pyplot as plt

def save(model,soln_index):
    """ Save new parameters """
    import analysis_scripts.plot_depsilon_native as eplot
    
    cwd = os.getcwd()
    eps0 = model.pairwise_strengths
    deps = np.loadtxt("xp_%d.dat" % soln_index)
    fitit = 1
    
    ##calculate the scaling based upon what would actually happen: if it's already at 0.01, and goes more negative, now it effectively shows a deps there of 0, and not some arbitrarily large number
    neweps_effective = eps0 + deps
    neweps_effective[neweps_effective < 0.01] = 0.01
    deps_effective = neweps_effective - eps0
    
    
    factor = np.linalg.norm(deps_effective)/np.linalg.norm(eps0)
    max_step = np.max(np.abs(deps_effective/eps0))
    max_step_factor = 0.3
        
    if factor > 0.3:
        deps = (deps*max_step_factor) / max_step
        print "Scaling down to 0.3 by maximum step"
        fitit = max_step_factor/max_step
    
    eplot.plot_epsilons_bin(deps,"d-epsilon",model)
    eplot.plot_epsilons(deps,"d-epsilon",model)
    
    neweps = eps0 + deps
    neweps[neweps < 0.01] = 0.01
    
    model.update_model_param_values(neweps)
    
    plt.figure()
    plt.hist(neweps, 50, alpha=0.75)
    plt.xlabel("epsilon",fontsize=20)
    plt.ylabel("number",fontsize=20)
    plt.title("spared of epsilons", fontsize=20)
    plt.savefig("eps_spread.png")
    
    #estimate_lambda()
    open("%s/pairwise_params" % cwd,"w").write(model.pairwise_param_file_string)
    open("%s/model_params" % cwd,"w").write(model.model_param_file_string)
    open("%s/fitting_scale" % cwd,"w").write("%f"%fitit)
    model.pairwise_params_file_location = "%s/pairwise_params" % cwd
    model.model_params_file_location = "%s/model_params" % cwd

def estimate_lambda():
    svs = np.loadtxt("singular_values.dat")
    index = 0
    num = np.shape(svs)[0]
    cwd = os.getcwd()
    for i in range(num-1):
        if svs[i]/svs[i+1] > 10000:
            index = num - i - 1
    open("%s/Lambda_index.txt"%cwd, "w").write("%d"%index)
    print "biggest difference is for a value of %f and %f" % (svs[num-1-index], svs[num-index]) 
        
        
        
        
        
        
        
        
        
        
        
    

