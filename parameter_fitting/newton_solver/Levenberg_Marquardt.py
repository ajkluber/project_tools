''' Solve for new parameters

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004. Solutions are
found with the Levenberg-Marquardt algorithm, an algorithm for damped 
least-squares.

References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
'''

import matplotlib.pyplot as plt
import numpy as np

import save_and_plot

def find_solutions(model, scaling=False, chosen_cutoffs=None, simplify=False):
    ''' Solve for new parameters using the Levenberg-Marquardt algorithm 

    Description:
    
        Solves the ill-posed (rectangular) system for the new parameters
    using the damped least-squares, Levenberg-Marquardt algorithm. Ill-posed
    systems usually have some small singular values that would blow the 
    solution up if the matrix inverse was calculated. The idea is to reduce 
    the influence of very small singular values, which have a lower signal 
    to noise ratio, on the solution. 

        The usual least-squares solution seeks to minimize the norm of the
    residual, 
                    Jx = df   =>  x s.t. min||Jx - df||

    damped least-squares modifies the problem to,

                A = J^T*J

    (A + lambda*diag(A))x = J^T*df  =>  x s.t. min{||Jx - b|| + lambda*||x||}

    Now a solution is sought that simultaneously minimizes the residual 
    and the norm of the solution. The choice of the damping parameter lambda
    can be determined

    '''

    target_feature = np.loadtxt("target_feature.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    J = np.loadtxt("Jacobian.dat")

    df = target_feature - sim_feature

    #if scaling:
    #    JTdf = np.dot(J.T,df)
    #    JTJ = np.dot(J.T,J)
    #    n_rows = JTJ.shape[0]
    #    ## Scale the diagonal by the curvature. This makes it the Levenberg-Marquardt method
    #    Levenberg = np.identity(n_rows)
    #    Levenberg[(np.arange(n_rows),np.arange(n_rows))] = np.diag(JTJ)

    #u,s,v = np.linalg.svd(JTJ)
    u,s,v = np.linalg.svd(J)
    if np.log10(min(s)) < -10:
        smin = -10
    else:
        smin = np.log10(min(s)) - 2
    ## The damping parameter, lambda, is scanned from below the smallest
    ## singular value to above the largest singular value.
    ## The contributions from singular values less than the damping 
    ## parameter are 'damped' out by a filter function.
    Lambdas = np.logspace(smin,np.log10(max(s)),num=200)

    norm_eps = np.linalg.norm(model.model_param_values[model.fitting_params])
    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    for i in range(len(Lambdas)):
        S = np.zeros(J.shape) 
        Lambda = Lambdas[i]

        #S[np.arange(len(s)),np.arange(len(s))] = s/((s**2) + (Lambda**2))
        S[np.diag_indices(len(s))] = s/((s**2) + (Lambda**2))

        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_soln = np.dot(J_pinv,df)

        residual = np.dot(J,x_soln) - df

        nrm_soln.append(np.linalg.norm(x_soln))
        nrm_resd.append(np.linalg.norm(residual))
        solutions.append(x_soln)

        #J_use = np.dot(v.T,np.dot(S.T,u.T))     ## This isn't right
        #cond_num = np.linalg.norm(J_use)*np.linalg.norm(J_pinv)
        #condition_number.append(np.linalg.norm(lhs)*np.linalg.norm(np.linalg.inv(lhs)))
        #condition_number.append(np.linalg.norm(lhs)*np.linalg.norm(np.linalg.inv(lhs)))

    condition_number = []
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
    #model.contact_epsilons = np.ones(model.n_contacts,float)

    cwd = os.getcwd()
    if not os.path.exists("%s/iteration_%d/test" % (name,iteration)):
        os.mkdir("%s/iteration_%d/test" % (name,iteration))
    os.chdir("%s/iteration_%d/test" % (name,iteration))

    find_solutions(model)

    os.chdir(cwd)



