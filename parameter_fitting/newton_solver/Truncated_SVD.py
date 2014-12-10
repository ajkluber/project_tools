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

def find_solutions(model,method):

    target_feature = np.loadtxt("target_feature.dat")
    target_feature_err = np.loadtxt("target_feature_err.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    sim_feature_err = np.loadtxt("sim_feature_err.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    Jacobian_err = np.loadtxt("Jacobian_err.dat")
    J = Jacobian

    epsilons = model.contact_epsilons
    norm_eps = np.linalg.norm(epsilons)

    ## Normalize the target step. 
    df = target_feature - sim_feature

    u,s,v = np.linalg.svd(J)
    temp  = list(0.5*(s[:-1] + s[1:])) + [0.0]
    temp.reverse()
    cutoffs = np.array(temp)

    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    Taus = []
    for i in range(len(cutoffs)):
        S = np.zeros(J.shape) 
        tau = cutoffs[i] 
        s_use = s[s > tau]
        S[np.arange(len(s_use)),np.arange(len(s_use))] = 1./s_use
        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_soln = np.dot(J_pinv,df)  ## particular

        #x_soln = np.dot(v.T,np.dot(S.T,np.dot(u.T,df)))
        if cplex:
            solution_cplex = func(J,cutoff,weight)

        residual = np.dot(J,x_soln) - df

        nrm_soln.append(np.linalg.norm(x_soln))
        nrm_resd.append(np.linalg.norm(residual))
        solutions.append(x_soln)
        Taus.append(tau)

        J_use = np.dot(v.T,np.dot(S.T,u.T))     ## This isn't right
        cond_num = np.linalg.norm(J_use)*np.linalg.norm(J_pinv)
        condition_number.append(cond_num)

    save_and_plot.save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s)


if __name__ == '__main__':
    import os
    import argparse
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    
    model = mdb.check_inputs.load_model("%s" % name,dry_run=True)
    model.iteration = iteration
    model.contact_epsilons = np.ones(model.n_contacts,float)
    method = "ddG_MC2004"

    cwd = os.getcwd()
    if not os.path.exists("%s/iteration_%d/test" % (name,iteration)):
        os.mkdir("%s/iteration_%d/test" % (name,iteration))
    os.chdir("%s/iteration_%d/test" % (name,iteration))

    find_solutions(model,method)

    os.chdir(cwd)



