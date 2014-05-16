import numpy as np
import os

import mdtraj as md
import cplex

import phi_values as phi

import model_builder.models as models
import model_builder.systems as systems

'''
April 30 2014
Alexander Kluber

Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004.
'''

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_thermodynamic_perturbation(Model,System,append_log,coord="Q"):
    ''' First task is to calculate the perturbations for each mutation for
        each frame in the trajectory.   May be generalized in the future or 
        moved inside Model to deal with Models with multiple parameters per
        interaction (e.g. desolvation barrier, etc.)
    '''
    
    #append_log(System.subdir,"Starting: Calculating_MC2004")

    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = phi.get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))

    if not os.path.exists(savedir+"/mut"):
        os.mkdir(savedir+"/mut")

    files = ["M.dat","ddG.dat","eps.dat"] 
    flag = np.array([ not os.path.exists(savedir+"/mut/"+file) for file in files ])
    if np.any(flag):
        print "  One of the following does not exist: M.dat, ddG.dat, eps.dat. Calculating."
        os.chdir(System.subdir)
        ddG, eps, M = calculate_matrix_ddG_eps_M(Model,System,savedir,beta,coord)
        os.chdir(cwd)
    else:
        print "  Loading M.dat, ddG.dat, eps.dat"
        ddG = np.loadtxt(savedir+"/mut/ddG.dat")
        eps = np.loadtxt(savedir+"/mut/eps.dat")
        M = np.loadtxt(savedir+"/mut/M.dat")

    for cutoff in [1.25]:
        LP_problem, solution, x_particular, N = apply_constraints_linear_objective(Model,System,savedir,ddG,eps,M,cutoff=cutoff)
        #LP_problem, solution, x_particular, N = apply_constraints_quadratic_objective(Model,System,savedir,ddG,eps,M,cutoff=cutoff)

        ## New Parameters
        epsilon_prime = eps + x_particular + np.dot(N,solution)

        np.savetxt(savedir+"/mut/line_delta_eps_"+str(cutoff)+".dat",x_particular+np.dot(N,solution))
        np.savetxt(savedir+"/mut/line_eps_prime_"+str(cutoff)+".dat",epsilon_prime)
        #np.savetxt(savedir+"/mut/quad_delta_eps_"+str(cutoff)+".dat",x_particular+np.dot(N,solution))
        #np.savetxt(savedir+"/mut/quad_eps_prime_"+str(cutoff)+".dat",epsilon_prime)
    

    return LP_problem, solution, x_particular, eps, N
    #append_log(System.subdir,"Finished: Calculating_MC2004")

def apply_constraints_quadratic_objective(Model,System,savedir,ddG,eps,M,cutoff=0.5):
    ''' Search the nullspace dimensions for a solution that minimizes 
        the change in stability and keeps the contacts attractive. Uses
        cplex linear programming package.

        status 4-29-2014 WORKS!
    '''

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
    Mpinv = np.linalg.pinv(M,rcond=cutoff)
    x_particular = np.dot(Mpinv,ddG)
    np.savetxt(savedir+"/mut/x_p.dat",x_particular)
    #print x_particular     ## DEBUGGING

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## S[:M.shape[1],:M.shape[1]] = np.diag(s)
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
    u,s,v = np.linalg.svd(M)
    rank = len(s)

    ## Nullspace basis vectors are the last n-r columns of the matrix v.T. As a check
    ## all entries of the matrix np.dot(M,N) should be extremely small ~0. Because 
    ## they've been sent to the nullspace.
    N = v.T[:,M.shape[0]:]

    print N.shape
    ## Objective coefficients are sum of nullspace vectors. This comes from
    ## requiring the same average contact strength.
    #objective_coeff = list(sum(N))
    objective_coeff = list(2.*np.dot(x_particular.T,N))

    ## Further constraints come from requiring that native contacts remain
    ## attractive.
    right_hand_side = list(eps + x_particular)
    column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    row_names = [ "c"+str(i) for i in range(N.shape[0]) ]
    rows = [ [column_names,list(-N[i,:])]  for i in range(len(N)) ]
    senses = "L"*len(right_hand_side)

    ## Set quadratic terms in objective.
    objective_quadratic_coefficients = [ 1. for j in range(N.shape[1]) ]

    ## DEBUGGING
    #print len(objective_coeff),len(names)
    #print right_hand_side
    #print rows
    #print rows[0]
    #print senses
    #print N
    #print len(rows), len(senses), len(rows[0][0]), len(rows[0][1])

    ## Set upper and lower bounds on the solution. Arbitrary. Hopefullly these 
    ## don't matter.
    upper_bounds = list(10000.*np.ones(N.shape[1]))
    lower_bounds = list(-10000.*np.ones(N.shape[1]))

    ## DEBUGGING
    #print len(upper_bounds)

    ## Populate cplex linear programming problem
    LP_problem = cplex.Cplex()
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
    LP_problem.objective.set_quadratic(objective_quadratic_coefficients)

    ## Let cplex do the hard work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution_lambda = LP_problem.solution.get_values()

    print "Cplex summary:"
    print "status: ",status
    print "solution:",solution_lambda
    
    #epsilon_new = eps + np.dot(N,np.array(solution))

    return LP_problem, solution_lambda, x_particular, N


def apply_constraints_linear_objective(Model,System,savedir,ddG,eps,M,cutoff=0.5):
    ''' Search the nullspace dimensions for a solution that minimizes 
        the change in stability and keeps the contacts attractive. Uses
        cplex linear programming package.

        status 4-29-2014 WORKS.
    '''

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
    Mpinv = np.linalg.pinv(M,rcond=cutoff)
    x_particular = np.dot(Mpinv,ddG)
    np.savetxt(savedir+"/mut/x_p.dat",x_particular)
    #print x_particular     ## DEBUGGING

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## S[:M.shape[1],:M.shape[1]] = np.diag(s)
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
    u,s,v = np.linalg.svd(M)
    rank = len(s)

    ## Nullspace basis vectors are the last n-r columns of the matrix v.T. As a check
    ## all entries of the matrix np.dot(M,N) should be extremely small ~0. Because 
    ## they've been sent to the nullspace.
    N = v.T[:,M.shape[0]:]

    print N.shape
    ## Objective coefficients are sum of nullspace vectors. This comes from
    ## requiring the same average contact strength.
    objective_coeff = list(sum(N))

    ## Further constraints come from requiring that native contacts remain
    ## attractive.
    right_hand_side = list(eps + x_particular)
    column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    row_names = [ "c"+str(i) for i in range(N.shape[0]) ]
    rows = [ [column_names,list(-N[i,:])]  for i in range(len(N)) ]
    senses = "L"*len(right_hand_side)

    ## DEBUGGING
    #print len(objective_coeff),len(names)
    #print right_hand_side
    #print rows
    #print rows[0]
    #print senses
    #print N
    #print len(rows), len(senses), len(rows[0][0]), len(rows[0][1])

    ## Set upper and lower bounds on the solution. Arbitrary. Hopefullly these 
    ## don't matter.
    upper_bounds = list(10000.*np.ones(N.shape[1]))
    lower_bounds = list(-10000.*np.ones(N.shape[1]))

    ## DEBUGGING
    #print len(upper_bounds)

    ## Populate cplex linear programming problem
    LP_problem = cplex.Cplex()
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)

    ## Let cplex do the hard work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution_lambda = LP_problem.solution.get_values()

    print "Cplex summary:"
    print "status: ",status
    print "solution:",solution_lambda
    
    #epsilon_new = eps + np.dot(N,np.array(solution))

    return LP_problem, solution_lambda, x_particular, N

    ## Save new parameters in a BeadBead.dat file and indicate to use the path
    ## to the file in the System object.
    
def calculate_matrix_ddG_eps_M(Model,System,savedir,beta,coord):
    ''' Calculates and saves and returns the matrix from equation (9) in 
        Matysiak Clementi 2004. '''

    print "  Getting state bounds for coordinate:",coord
    bounds, states = phi.get_state_bounds(savedir,coord)
    num_states = len(states)

    print "  Loading mutants"
    mutants = [ x.split()[1]+x.split()[0]+x.split()[2] for x in open("mutants/mutations.dat","r").readlines()[1:] ]

    print "  Loading trajectory, epsij, deltaij, sigij"
    sigij,epsij,deltaij,interaction_nums,keep_interactions,pairs,traj,traj_dist = phi.load_eps_delta_sig_traj(savedir)
    Fij = phi.get_mutant_fij(mutants,keep_interactions)
    qij = phi.get_Qij(Model,traj_dist,sigij,deltaij,interaction_nums)
    print "  Loading dH for mutants"
    dH = phi.get_mutant_dH(savedir,mutants)

    np.savetxt(savedir+"/mut/eps.dat",epsij)
    np.savetxt(savedir+"/mut/delta.dat",deltaij)
    np.savetxt(savedir+"/mut/sigma.dat",sigij)

    ## DEBUGGING 
    #print qij.shape, dH.shape, Fij.shape
    #print qij.shape, dH.shape
    #print "Qij shape",qij[states[1],:].shape
    #print "dH shape",dH[:,states[1]].shape
    #print "dH shape",sum(dH[:,states[1]])

    ## Load ddG from experiment and theory. 
    print "  Loading ddG from simulation"
    ddGsim_TS_D, ddGsim_N_D = phi.get_sim_ddG(savedir,coord,bounds)
    print "  Loading ddG from experiment"
    ddGexp_TS_D, ddGexp_N_D = phi.get_exp_ddG()
    ddG_all = np.concatenate(((ddGexp_TS_D - ddGsim_TS_D),(ddGexp_N_D - ddGsim_N_D)), axis=0)

    np.savetxt(savedir+"/mut/ddG.dat",ddG_all)
    
    ## Then compute each term of MC2004 equation (9) in a vectorized fashion. 
    ## This is the formula derived from thermodynamic perturbation.
    termA_TS_D = np.ones(Fij.shape,float)*(sum(qij[states[1],:])/float(len(qij[states[1],:])) - sum(qij[states[0],:])/float(len(qij[states[0],:])))
    termA_N_D = np.ones(Fij.shape,float)*(sum(qij[states[2],:])/float(len(qij[states[2],:])) - sum(qij[states[0],:])/float(len(qij[states[0],:])))
    termA_all = np.concatenate((termA_TS_D,termA_N_D),axis=0)
    
    termB_all = 1. - np.concatenate((Fij,Fij),axis=0) 

    termC_TS_D = ((np.dot(np.exp(beta*dH[:,states[1]]),qij[states[1],:])).T/sum(np.exp(beta*dH[:,states[1]]).T)).T \
              - ((np.dot(np.exp(beta*dH[:,states[0]]),qij[states[0],:])).T/sum(np.exp(beta*dH[:,states[0]]).T)).T
    termC_N_D = ((np.dot(np.exp(beta*dH[:,states[2]]),qij[states[2],:])).T/sum(np.exp(beta*dH[:,states[2]]).T)).T \
              - ((np.dot(np.exp(beta*dH[:,states[0]]),qij[states[0],:])).T/sum(np.exp(beta*dH[:,states[0]]).T)).T
    termC_all = np.concatenate((termC_TS_D,termC_N_D),axis=0)

    ## Then contruct the matrix M. 
    M = -beta*( termA_all -  termB_all*termC_all)

    ## DEBUGGING
    #print "Term A", termA_all.shape
    #print "Term B", termB_all.shape
    #print "Term C", termC_all.shape
    #print "M", M.shape
    #print "ddG", ddG_all.shape

    np.savetxt(savedir+"/mut/termA.dat",termA_all)
    np.savetxt(savedir+"/mut/termB.dat",termB_all)
    np.savetxt(savedir+"/mut/termC.dat",termC_all)
    np.savetxt(savedir+"/mut/M.dat",M)
    
    return ddG_all,epsij,M

if __name__ == '__main__':
    ## TESTING calculating ddG, phi values.
    def dummy_func(sub,string):
        pass 
    
    subdirs = ["r15"]
    Models = models.load_models(subdirs,dryrun=True)
    Systems = systems.load_systems(subdirs)
    Model = Models[0]
    System = Systems[0]
    path = System.subdir+"/"+System.mutation_active_directory+"/131.17_agg"

    '''
    #bounds, states = phi.get_state_bounds(path,"Q") ## DEBUGGING
    dH, states = calculate_phi_values(Model,System,dummy_func)
    '''

    LP_problem, solution, x_particular, eps, N = calculate_thermodynamic_perturbation(Model,System,dummy_func)
