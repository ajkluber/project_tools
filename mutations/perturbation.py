""" Apply thermodynamic perturbation via MC2004 algorithm

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import numpy as np
import os

import mdtraj as md
import cplex

import phi_values as phi
import mutatepdbs as mut

import model_builder.models as models
import model_builder.systems as systems

import matplotlib.pyplot as plt

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_MC2004_perturbation(Model,System,append_log,coord="Q",newbeadbead="NewBeadBead.dat",target_ratio=0.95):
    """ Calculate the new contact parameters with Matysiak Clementi 2004 method

    Description:

        Use Matysiak, Clementi 2004 perturbation technique to solve for new
    contact parameters. See reference (1) for more details. The procedure 
    is based off of trying to match experimental DeltaDelta G's (from 
    phi value analysis) to simulation DeltaDelta G's. It involves Taylor 
    expanding the simulation DeltaDelta G's around 


    Reference:

    (1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
    the Characterization of the Protein Folding Landscape of S6: How Far Can a
    Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
    """
    
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
        print "  IMPORTANT: Look at the singular value spectrum and choose a number of singular values to use."
        print "  IMPORTANT: Make sure ",savedir+"/mut/num_singular_values_include.txt  exists."
        os.chdir(cwd)
    else:
        print "  Loading M.dat, ddG.dat, eps.dat"
        ddG = np.loadtxt(savedir+"/mut/ddG.dat")
        eps = np.loadtxt(savedir+"/mut/eps.dat")
        M = np.loadtxt(savedir+"/mut/M.dat")

    u,s,v = np.linalg.svd(M)
    s_norm = s/max(s)
    cutoffs = s_norm - 0.01*np.ones(len(s_norm))

    if not os.path.exists(savedir+"/mut/num_singular_values_include.txt"):
        print "ERROR!"
        print "  Need file ", savedir+"/mut/num_singular_values_include.txt", " to proceed"
        print "  Exiting"
        raise SystemExit
    else:
        temp = open(savedir+"/mut/num_singular_values_include.txt").read()[:-1]
        if temp.endswith("xp"):
            compute_xp = True
            num_singular_values = int(temp.split()[0])
            savebeadbeadxp = savedir+"/mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = savebeadbeadxp
        else:
            compute_xp = False
            num_singular_values = int(temp)
            savebeadbeadxp = savedir+"/mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = savedir+"/mut/"+newbeadbead
        cutoff = cutoffs[num_singular_values]
        print "  Using ",num_singular_values," singular values. Cutoff of = ",cutoff

    append_log(System.subdir,"Starting: Calculating_MC2004") 

    if compute_xp == True:
        print "  Using ONLY the particular solution x_p as delta_eps!"
        Mpinv = np.linalg.pinv(M,rcond=cutoff)
        x_particular = np.dot(Mpinv,ddG)
        np.savetxt(savedir+"/mut/x_p.dat",x_particular)

        delta_eps = x_particular
        ratio = np.linalg.norm(delta_eps)/np.linalg.norm(eps)

        delta_eps_xp = x_particular
        ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)
    else:
        print "  Applying CPLEX to the particular solution to get delta_eps!"
        LP_problem, solution, x_particular, N = apply_constraints_with_cplex(Model,System,savedir,ddG,eps,M,cutoff)

        print "    Solution found!"
        delta_eps = x_particular + np.dot(N,solution)
        ratio = np.linalg.norm(delta_eps)/np.linalg.norm(eps)

        delta_eps_xp = x_particular
        ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)

    np.savetxt(savedir+"/mut/delta_eps.dat",delta_eps)
    np.savetxt(savedir+"/mut/delta_eps_xp.dat",delta_eps_xp)
    print "    Norm of perturbation, |deps|/|eps| = ", ratio
    print "    Norm of just x_p perturbation, |x_p|/|eps| = ", ratio_xp

    print "  Saving new parameters as: ",savebeadbead
    save_new_parameters(sub,eps,delta_eps,delta_eps_xp,savebeadbead, savebeadbeadxp, T)

    Model.contact_energies = savebeadbead

    append_log(System.subdir,"Finished: Calculating_MC2004") 

def apply_constraints_with_cplex(Model,System,savedir,ddG,eps,M,cutoff):
    """ Construct and solve a linear/quadratic programming problem for new parameters.

    Description:

        Use the IMB ILOG Cplex optimization library to search the nullspace 
    dimensions for a solution that satisfies the given constraints and mimimizes
    the given objective function

    """

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
    print "  Matrix M singular values saved as: ", savedir+"/mut/singular_values.dat (and as ..._norm.dat)"
    print "  Normed singular value spectrum:"
    print s/max(s)
    np.savetxt(savedir+"/mut/singular_values.dat",s)
    np.savetxt(savedir+"/mut/singular_values_norm.dat",s/max(s))

    ## Nullspace basis vectors are the last n-r columns of the matrix v.T. As a check
    ## all entries of the matrix np.dot(M,N) should be extremely small ~0. Because 
    ## they've been sent to the nullspace.
    N = v.T[:,M.shape[0]:]

    ############# OBJECTIVE FUNCTION ###############
    ## Objective coefficients are sum of nullspace vectors. This comes from
    ## requiring the same average contact strength.

    ## Objective 1: Linear objective that seeks to minimize the change in 
    ## average contact strength. Found to give crazy results.
    #objective_coeff = list(sum(N))

    ## Objective 2: Quadratic objective that seeks to minimize the size of the
    ## perturbation.
    objective_coeff = list(2.*np.dot(x_particular.T,N))

    ############# CONSTRAINTS ###############
    ## Linear constraints are applied to keep the resulting parameters within
    ## a certain range. Only uncomment one version of constraints at a time.
    
    ## Constraints version 1:
    ## Require that native contacts remain attractive, 
    ## i.e. eps'_ij in the interval (0,inf)
    #right_hand_side = list(eps + x_particular)
    #column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    #row_names = [ "c"+str(i) for i in range(N.shape[0]) ]
    #rows = [ [column_names,list(-N[i,:])]  for i in range(len(N)) ]
    #senses = "L"*len(right_hand_side)

    ## Constraints version 2:
    ## Require that native contacts remain attractive, with lower bound
    ## i.e. eps'_ij in the interval (lower_bound,inf)
    #a = 0.10
    #right_hand_side = list(eps + x_particular - a*np.ones(len(eps)) )
    #column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    #row_names = [ "c"+str(i) for i in range(N.shape[0]) ]
    #rows = [ [column_names,list(-N[i,:])]  for i in range(len(N)) ]
    #senses = "L"*len(right_hand_side)

    ## Constraints version 3:
    ## Require that native contacts remain within a range of their previous
    ## values
    ## i.e. eps'_ij in the interval (eps_ij/a, a*eps_ij) for some number a.
    #a = 3.
    #right_hand_side = list((1./a)*eps + x_particular)
    #right_hand_side_2 = list((1.-a)*eps + x_particular)
    #right_hand_side.extend(right_hand_side_2)
    #column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    #row_names = [ "c"+str(i) for i in range(2*N.shape[0]) ]
    #senses = "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)
    #rows = []
    #for n in range(2):
    #    for i in range(len(N)):
    #        rows.append([ column_names, list(-N[i,:]) ])

    ## Constraints version 4:
    ## Require that native contacts remain within a fixed interval.
    eps_lower_bound = 0.1
    eps_upper_bound = 4.0
    right_hand_side = list(eps + x_particular - eps_lower_bound*np.ones(len(eps)))
    right_hand_side_2 = list(eps + x_particular - eps_upper_bound*np.ones(len(eps)))
    right_hand_side.extend(right_hand_side_2)

    column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    row_names = [ "c"+str(i) for i in range(2*N.shape[0]) ]
    senses = "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)

    rows = []
    for n in range(2):
        for i in range(len(N)):
            rows.append([ column_names, list(-N[i,:]) ])

    ## Set quadratic terms in objective.
    objective_quadratic_coefficients = [ 1. for j in range(N.shape[1]) ]

    ## Set upper and lower bounds on the solution. Arbitrary. Hopefullly these 
    ## don't matter. These are bounds on vector lambda
    upper_bounds = list(10000.*np.ones(N.shape[1]))
    lower_bounds = list(-10000.*np.ones(N.shape[1]))

    ## Populate cplex linear programming problem
    LP_problem = cplex.Cplex()
    LP_problem.set_results_stream(None)
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
    LP_problem.objective.set_quadratic(objective_quadratic_coefficients)

    ## Let cplex do the hard work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution_lambda = LP_problem.solution.get_values()

    ## Print cplex summary
    #print "Cplex summary:"
    #print "status: ",status
    #print "solution:",solution_lambda
    
    return LP_problem, solution_lambda, x_particular, N

def calculate_matrix_ddG_eps_M(Model,System,savedir,beta,coord):
    ''' Calculates and saves and returns the matrix from equation (9) in 
        Matysiak Clementi 2004. '''

    print "  Getting state bounds for coordinate:",coord
    bounds, states = phi.get_state_bounds(savedir,coord)
    num_states = len(states)

    print "  Loading mutants"
    os.chdir("mutants")
    mut_indx,wt_res,mut_res = mut.get_core_mutations()
    print "  Loading ddG from experiment"
    ddGexp_N_D,ddGexp_N_D_err,ddGexp_TS_D,ddGexp_TS_D_err = mut.get_core_mutation_ddG()
    os.chdir("..")
    mutants = [ wt_res[i]+mut_indx[i]+mut_res[i] for i in range(mut_indx.shape[0]) ]
    #print "  Computing perturbation with the following mutants"
    #print mutants

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

    ## Load ddG from theory. 
    print "  Loading ddG from simulation"
    ddGsim_TS_D, ddGsim_N_D = mut.get_sim_ddG(savedir,mutants,coord,bounds)

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

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## S[:M.shape[1],:M.shape[1]] = np.diag(s)
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
    u,s,v = np.linalg.svd(M)
    rank = len(s)
    np.savetxt(savedir+"/mut/singular_values.dat",s)
    np.savetxt(savedir+"/mut/singular_values_norm.dat",s/max(s))
    print "  Matrix M singular values saved as: ", savedir+"/mut/singular_values.dat"
    print "  Normed singular value spectrum:"
    print s/max(s)

    np.savetxt(savedir+"/mut/termA.dat",termA_all)
    np.savetxt(savedir+"/mut/termB.dat",termB_all)
    np.savetxt(savedir+"/mut/termC.dat",termC_all)
    np.savetxt(savedir+"/mut/M.dat",M)
    print "  Matrices: termA, termB, termC, and M computed and saved to ", savedir+"/mut"
    
    return ddG_all,epsij,M

def save_new_parameters(sub,eps,delta_eps,delta_eps_xp,savebeadbead,savebeadbeadxp, T):
    """ Save new parameters as a BeadBead file

    Description:

    """
    beadbead, keep_interactions = phi.load_beadbead(sub+"/"+T+"_1")

    ## Saving new parameters
    epsilon_prime = eps + delta_eps
    epsij = beadbead[:,6].astype(float)
    epsij[keep_interactions != 0] = epsilon_prime
    
    beadbead_string = ''
    for rownum in range(beadbead.shape[0]): 
        i_idx = int(beadbead[rownum,0])
        j_idx = int(beadbead[rownum,1])
        resi_id = beadbead[rownum,2]
        resj_id = beadbead[rownum,3]
        interaction_num = str(beadbead[rownum,4])
        sig = float(beadbead[rownum,5])
        Knb = epsij[rownum]
        delta = float(beadbead[rownum,7])
        beadbead_string += '%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n' % \
                (i_idx,j_idx,resi_id,resj_id,interaction_num,sig,Knb,delta)

    open(savebeadbead,"w").write(beadbead_string)

    ## Also saving the x_p solution as well for good measure.
    epsilon_prime = eps + delta_eps_xp
    epsij = beadbead[:,6].astype(float)
    epsij[keep_interactions != 0] = epsilon_prime

    beadbead_string = ''
    for rownum in range(beadbead.shape[0]): 
        i_idx = int(beadbead[rownum,0])
        j_idx = int(beadbead[rownum,1])
        resi_id = beadbead[rownum,2]
        resj_id = beadbead[rownum,3]
        interaction_num = str(beadbead[rownum,4])
        sig = float(beadbead[rownum,5])
        Knb = epsij[rownum]
        delta = float(beadbead[rownum,7])
        beadbead_string += '%5d%5d%8s%8s%5s%16.8E%16.8E%16.8E\n' % \
                (i_idx,j_idx,resi_id,resj_id,interaction_num,sig,Knb,delta)

    open(savebeadbeadxp,"w").write(beadbead_string)

if __name__ == '__main__':
    def dummy_func(sub,string):
        pass 
    
    subdirs = ["r16"]
    Tf_choice = open(subdirs[0]+"/Mut_0/Tf_choice.txt","r").read()[:-1]
    Models = models.load_models(subdirs,dryrun=True)
    Systems = systems.load_systems(subdirs)
    Model = Models[0]
    System = Systems[0]
    path = System.subdir+"/"+System.mutation_active_directory+"/"+Tf_choice+"_agg"


    '''
    #bounds, states = phi.get_state_bounds(path,"Q") ## DEBUGGING
    dH, states = calculate_phi_values(Model,System,dummy_func)
    '''
    target_ratio = 1.0
    calculate_MC2004_perturbation_debug(Model,System,dummy_func,newbeadbead="test_threshold_sing_val_0.2_"+str(target_ratio)+".dat",target_ratio=target_ratio)

    raise SystemExit
    newbeadbead="test_threshold_lb_0.01_ub_4.0_"+str(target_ratio)+".dat"

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

    u,s,v = np.linalg.svd(M)
    s_norm = s/max(s)

    cutoffs = s_norm - 0.01*np.ones(len(s_norm))

    ## Linear search for best singular value cutoff between 0.0 and 0.5 
    ## Stopping criteria is when the perturbation is 'small' i.e. the
    ## ratio of the perturbation to the norm of the current parameters is
    ## less than 1, so around 0.95 ==> |depsilon|/|epsilon| ~ 0.95 < 1 

    tolerance = 0.04
    #target_ratio = 0.95
    ratio = 0
    iteration = 1 

    print "Iteration  Cutoff    Ratio "
    cutoff = 0.8
    xps = []

    norm_delta_eps = []
    std_delta_eps = []

    std_eps_prime = []
    avg_eps_prime = []
    for cutoff in cutoffs[1:len(cutoffs)/2]:
        LP_problem, solution, x_particular, N = apply_constraints_with_cplex(Model,System,savedir,ddG,eps,M,cutoff)
        delta_eps = x_particular + np.dot(N,solution)
        eps_prime = delta_eps + eps

        xps.append(x_particular)

        ratio = np.linalg.norm(delta_eps)/np.linalg.norm(eps)

        norm_delta_eps.append(ratio)

        std_eps_prime.append(np.std(eps_prime))
        avg_eps_prime.append(np.mean(eps_prime))

        print iteration, cutoff, ratio
        iteration += 1 




