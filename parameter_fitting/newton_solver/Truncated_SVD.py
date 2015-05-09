''' Solve for new parameters

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004. Solutions are
found with the Truncated Singular Value Decomposition (TSVD) algorithm

References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
'''

import matplotlib.pyplot as plt
import numpy as np

import save_and_plot

def find_solutions(model):
    target_feature = np.loadtxt("target_feature.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    J = Jacobian

    norm_eps = np.linalg.norm(model.model_param_values[model.fitting_params])
    df = target_feature - sim_feature

    u,s,v = np.linalg.svd(J)
    temp  = list(0.5*(s[:-1] + s[1:])) + [0.0]
    temp.reverse()
    Lambdas = np.array(temp)
    print s

    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    for i in range(len(Lambdas)):
        S = np.zeros(J.shape) 
        Lambda = Lambdas[i] 
        s_use = s[s > Lambda]
        #S[np.arange(len(s_use)),np.arange(len(s_use))] = 1./s_use
        S[np.diag_indices(len(s_use))] = 1./s_use
        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_soln = np.dot(J_pinv,df)  ## particular

        #x_soln = np.dot(v.T,np.dot(S.T,np.dot(u.T,df)))

        residual = np.dot(J,x_soln) - df

        nrm_soln.append(np.linalg.norm(x_soln))
        nrm_resd.append(np.linalg.norm(residual))
        solutions.append(x_soln)

        J_use = np.dot(v.T,np.dot(S.T,u.T))     ## This isn't right
        cond_num = np.linalg.norm(J_use)*np.linalg.norm(J_pinv)
        condition_number.append(cond_num)

    save_and_plot.save_solution_data(solutions,Lambdas,nrm_soln,nrm_resd,norm_eps,condition_number,s)


if __name__ == '__main__':
    import os
    import model_builder as mdb

    name = "S6" 
    model, fitopts = mdb.inputs.load_model("%s" % name)

    os.chdir("%s/iteration_%d/newton" % (name,0))
    find_solutions(model)
    os.chdir("../../..")
