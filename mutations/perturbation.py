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
import mutatepdbs as mute

import model_builder.models as models
import model_builder.systems as systems

import matplotlib.pyplot as plt

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_MC2004_perturbation(model,append_log,coord="Q",newbeadbead="NewBeadBead.dat",target_ratio=0.95):
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
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)
    T = phi.get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"

    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))

    if not os.path.exists(savedir+"/mut"):
        os.mkdir(savedir+"/mut")

    files = ["M.dat","ddG.dat","eps.dat"] 
    flag = np.array([ not os.path.exists(savedir+"/mut/"+file) for file in files ])
    if np.any(flag):
        print "  One of the following does not exist: M.dat, ddG.dat, eps.dat. Calculating."
        os.chdir(model.subdir)
        ddG, eps, M = calculate_matrix_ddG_eps_M(model,coord)
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
    #If s_norm is lower than 0.01, cutoffs should still be positve
    if cutoffs <0.:
        cutoffs = 0.

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
        cutoff = cutoffs[num_singular_values-1]
        print "  Using ",num_singular_values," singular values. Cutoff of = ",cutoff

    append_log(model.subdir,"Starting: Calculating_MC2004") 

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
        LP_problem, solution, x_particular, N = apply_constraints_with_cplex(model,savedir,ddG,eps,M,cutoff)

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

    model.contact_energies = savebeadbead

    append_log(model.subdir,"Finished: Calculating_MC2004") 

def apply_constraints_with_cplex(model,savedir,ddG,eps,M,cutoff):
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
    eps_lower_bound = 0.001
    eps_upper_bound = 20.0
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

def calculate_matrix_ddG_eps_M(model,coord):
    ''' Calculates and saves and returns the matrix from equation (9) in 
        Matysiak Clementi 2004. '''

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)

    os.chdir(model.subdir+"/mutants")
    print "  Loading mutants"
    mutants = mute.get_core_mutations()
    Fij, Fij_pairs, Fij_conts = phi.get_mutant_fij(model,mutants)
    print "  Loading ddG from experiment"
    ddGexp_N_D,ddGexp_N_D_err,ddGexp_TS_D,ddGexp_TS_D_err = mute.get_core_mutation_ddG()

    os.chdir(sub)

    bounds, states = phi.get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    ## 
    #print "  Loading ddG from simulation"
    #ddGsim_TS_D, ddGsim_N_D = mut.get_sim_ddG(mutants,coord)

    epsilons = model.contact_epsilons
    deltas = model.contact_deltas
    sigmas = model.contact_sigmas

    ## Loop over temperature directories.
    for n in range(len(directories)):
        T = temperatures[n]
        Tdir = directories[n]
        print " Computing matrix M for ", Tdir
        os.chdir(Tdir)
        if not os.path.exists("mut"):
            os.mkdir("mut")

        ## Boolean arrays that indicate which state each frame is in.
        ## States are defined by their boundaries along coordinate Q.
        Q = np.loadtxt("Q.dat")
        U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
        TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
        N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
        Nframes  = float(sum(N.astype(int)))
        Uframes  = float(sum(U.astype(int)))
        TSframes = float(sum(TS.astype(int)))

        ## Compute M vectorized
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        traj = md.load("traj.xtc",top="Native.pdb")
        rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape))

        ## Only count values of potential energy function where interaction is
        ## attractive.
        x = sigmas/rij
        x[(x > 1.09)] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero.
        Vij = epsilons*(5.*(x**12) - 6.*deltas*(x**10))
    
        Vij_U  = sum(Vij[U,:])/Uframes
        Vij_TS = sum(Vij[TS,:])/TSframes
        Vij_N  = sum(Vij[N,:])/Nframes

        #rij = md.compute_distances(traj,Fij_pairs[k])
        M = np.zeros((2*len(mutants),model.n_contacts),float)

        ## Compute the rows of matrix in equation (9) of reference (1)
        for k in range(len(mutants)):
            mut = mutants[k]
            print "   computing row for mutant ", mut
            ## Load dH_k for mutation. Eq. (2) minus Eq. (1) from reference (1).
            dHk = np.loadtxt("dH_"+mut+".dat")

            ## Thermal averages for matrix equation (9).
            expdHk_U  = sum(np.exp(-beta*dHk[U]))/Uframes
            expdHk_TS = sum(np.exp(-beta*dHk[TS]))/TSframes
            expdHk_N  = sum(np.exp(-beta*dHk[N]))/Nframes

            Vij_expdHk_U  = sum((Vij[U,:].T*np.exp(-beta*dHk[U])).T)/Uframes
            Vij_expdHk_TS = sum((Vij[TS,:].T*np.exp(-beta*dHk[TS])).T)/TSframes
            Vij_expdHk_N  = sum((Vij[N,:].T*np.exp(-beta*dHk[N])).T)/Nframes

            ## Compute all columns with Fij_k zero.
            M[k,:] = -beta*((Vij_TS - Vij_U) -  ((Vij_expdHk_TS/expdHk_TS) - (Vij_expdHk_U/expdHk_U)))
            M[k + len(mutants),:] = -beta*((Vij_N - Vij_U)  -  ((Vij_expdHk_N/expdHk_N) - (Vij_expdHk_U/expdHk_U)))

            ## Replace columns for which Fij_k is not zero.
            M[k,Fij_conts[k]] = -beta*((Vij_TS[Fij_conts[k]] - Vij_U[Fij_conts[k]])  - \
                (1. - Fij[k])*((Vij_expdHk_TS[Fij_conts[k]]/expdHk_TS) - \
                               (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))

            M[k + len(mutants),Fij_conts[k]] = -beta*((Vij_N[Fij_conts[k]] - Vij_U[Fij_conts[k]])  -  \
                (1. - Fij[k])*((Vij_expdHk_N[Fij_conts[k]]/expdHk_N)   - \
                              (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))
                
        u,s,v = np.linalg.svd(M)
        rank = len(s)
        np.savetxt("mut/M.dat",M)
        np.savetxt("mut/singular_values.dat",s)
        np.savetxt("mut/singular_values_norm.dat",s/max(s))
        plt.figure()
        plt.plot(s/max(s),'ro')
        plt.title(model.subdir+" "+Tdir)
        plt.savefig("mut/spectrum.pdf")
        plt.close()

        os.chdir("..")

    os.chdir(cwd)

    #return ddG_all,epsij,M

def save_new_parameters(sub,eps,delta_eps,delta_eps_xp,savebeadbead,savebeadbeadxp, T):
    """ Save new parameters as a BeadBead file

    Description:

    """
    beadbead, keep_interactions = phi.load_beadbead(sub+"/"+T+"_agg")

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
    
    subdirs = ["r15"]
    Models = models.load_models(subdirs,dryrun=True)
    model = Models[0]

    calculate_matrix_ddG_eps_M(model,"Q")
    raise SystemExit
    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)

    os.chdir(model.subdir+"/mutants")
    print "  Loading mutants"
    mutants = mute.get_core_mutations()
    Fij, Fij_pairs, Fij_conts = phi.get_mutant_fij(model,mutants)
    print "  Loading ddG from experiment"
    ddGexp_N_D,ddGexp_N_D_err,ddGexp_TS_D,ddGexp_TS_D_err = mute.get_core_mutation_ddG()

    os.chdir(sub)

    bounds, states = phi.get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    ## 
    #print "  Loading ddG from simulation"
    #ddGsim_TS_D, ddGsim_N_D = mut.get_sim_ddG(mutants,coord)

    epsilons = model.contact_epsilons
    deltas = model.contact_deltas
    sigmas = model.contact_sigmas

    ## Loop over temperature directories.
    for n in range(len(directories)):
        T = temperatures[n]
        Tdir = directories[n]
        print " Computing matrix M for ", Tdir
        os.chdir(Tdir)
        if not os.path.exists("mut"):
            os.mkdir("mut")

        ## Boolean arrays that indicate which state each frame is in.
        ## States are defined by their boundaries along coordinate Q.
        Q = np.loadtxt("Q.dat")
        U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
        TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
        N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
        Nframes  = float(sum(N.astype(int)))
        Uframes  = float(sum(U.astype(int)))
        TSframes = float(sum(TS.astype(int)))

        ## Compute M vectorized
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        traj = md.load("traj.xtc",top="Native.pdb")
        rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape))

        ## Only count values of potential energy function where interaction is
        ## attractive.
        x = sigmas/rij
        x[(x > 1.09)] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero.
        Vij = epsilons*(5.*(x**12) - 6.*deltas*(x**10))
    
        Vij_U  = sum(Vij[U,:])/Uframes
        Vij_TS = sum(Vij[TS,:])/TSframes
        Vij_N  = sum(Vij[N,:])/Nframes

        #rij = md.compute_distances(traj,Fij_pairs[k])
        M = np.zeros((2*len(mutants),model.n_contacts),float)

        ## Compute the rows of matrix in equation (9) of reference (1)
        for k in range(len(mutants)):
            mut = mutants[k]
            print "   computing row for mutant ", mut
            ## Load dH_k for mutation. Eq. (2) minus Eq. (1) from reference (1).
            dHk = np.loadtxt("dH_"+mut+".dat")

            ## Thermal averages for matrix equation (9).
            expdHk_U  = sum(np.exp(-beta*dHk[U]))/Uframes
            expdHk_TS = sum(np.exp(-beta*dHk[TS]))/TSframes
            expdHk_N  = sum(np.exp(-beta*dHk[N]))/Nframes

            Vij_expdHk_U  = sum((Vij[U,:].T*np.exp(-beta*dHk[U])).T)/Uframes
            Vij_expdHk_TS = sum((Vij[TS,:].T*np.exp(-beta*dHk[TS])).T)/TSframes
            Vij_expdHk_N  = sum((Vij[N,:].T*np.exp(-beta*dHk[N])).T)/Nframes

            ## Compute all columns with Fij_k zero.
            M[k,:] = -beta*((Vij_TS - Vij_U) -  ((Vij_expdHk_TS/expdHk_TS) - (Vij_expdHk_U/expdHk_U)))
            M[k + len(mutants),:] = -beta*((Vij_N - Vij_U)  -  ((Vij_expdHk_N/expdHk_N) - (Vij_expdHk_U/expdHk_U)))

            ## Replace columns for which Fij_k is not zero.
            M[k,Fij_conts[k]] = -beta*((Vij_TS[Fij_conts[k]] - Vij_U[Fij_conts[k]])  - \
                (1. - Fij[k])*((Vij_expdHk_TS[Fij_conts[k]]/expdHk_TS) - \
                               (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))

            M[k + len(mutants),Fij_conts[k]] = -beta*((Vij_N[Fij_conts[k]] - Vij_U[Fij_conts[k]])  -  \
                (1. - Fij[k])*((Vij_expdHk_N[Fij_conts[k]]/expdHk_N)   - \
                              (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))
                
        u,s,v = np.linalg.svd(M)
        rank = len(s)
        np.savetxt("mut/M.dat",M)
        np.savetxt("mut/singular_values.dat",s)
        np.savetxt("mut/singular_values_norm.dat",s/max(s))
        plt.figure()
        plt.plot(s/max(s),'ro')
        plt.title(model.subdir+" "+Tdir)
        plt.savefig("mut/spectrum.pdf")
        plt.close()

        os.chdir("..")

    os.chdir(cwd)

