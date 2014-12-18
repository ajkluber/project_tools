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
    
    open("%s/pairwise_params" % cwd,"w").write(model.pairwise_param_file_string)
    open("%s/model_params" % cwd,"w").write(model.model_param_file_string)
    model.contact_params_file_location = "%s/pairwise_params" % cwd
    model.model_params_file_location = "%s/model_params" % cwd

