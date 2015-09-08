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

import math

def find_solutions(model, scaling=False, chosen_cutoffs=None, simplify=False):
    cplex = False
    target_feature = np.loadtxt("target_feature.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    J = Jacobian

    epsilons = model.model_param_values
    norm_eps = np.linalg.norm(epsilons)

    ## Normalize the target step. 
    df = target_feature - sim_feature

    u,s,v = np.linalg.svd(J)
    
    #pick the lambdas (set of cutoffs)
    if chosen_cutoffs == None:
        temp  = list(0.5*(s[:-1] + s[1:])) + [0.0]
        temp.reverse()
        Lambdas = np.array(temp)
        
        #check and see length of Lambda, if it's super large, cut down by order of magnitude...:
        if len(Lambdas) > 200 and simplify:
            print "Number of cutoff values > 200, simplifying down to order of magnitudes"
            max_value = np.max(Lambdas)
            min_value = np.min(Lambdas[1:])
            max_power = math.floor(math.log(max_value, 10))
            min_power = math.floor(math.log(min_value, 10))
            Lambdas = np.zeros(2*(max_power - min_power)-1)
            for i in range(int(max_power - min_power)-1):
                Lambdas[2*i+1] = 10 ** (min_power + i +1)
                Lambdas[2*i+2] = 5*Lambdas[2*i+1]
            print Lambdas
    else:
        Lambdas = chosen_cutoffs
        
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
        if cplex:
            solution_cplex = func(J,cutoff,weight)

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



