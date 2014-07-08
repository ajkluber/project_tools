""" Apply thermodynamic perturbation via MC2004 algorithm

Description:

    Submodule to solve for the new parameter set using the thermodynamic
perturbation technique outlined in Matysiak Clementi 2004.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import matplotlib.pyplot as plt
import numpy as np
import os

import mdtraj as md
import cplex

import phi_values as phi
import mutatepdbs as mut

import model_builder.models as models
import model_builder.systems as systems

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_MC2004_perturbation(model,append_log,coord="Q",newbeadbead="NewBeadBead.dat"):
    """ Calculate the new contact parameters with Matysiak Clementi 2004 method

    Description:

        Use Matysiak, Clementi 2004 perturbation technique to solve for new
    contact parameters. See reference (1) for more details. The procedure 
    is based off of trying to match experimental DeltaDelta G's (from 
    phi value analysis) to simulation DeltaDelta G's. It involves Taylor 
    expanding the simulation DeltaDelta G's around 
        A linear system is solved for the parameter corrections and optionally
    adjusted with the nullspace degrees of freedom using a linear (quadratic)
    programming algorithm that seeks to minimize the norm of the resulting
    perturbation.

    Reference:

    (1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
    the Characterization of the Protein Folding Landscape of S6: How Far Can a
    Minimalist Model Go? J. Mol. Biol. 2004, 343, 235-248.
    """
    
    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)

    os.chdir(model.subdir+"/mutants")
    ddGexp, ddGexp_err = mut.get_exp_ddG()
    os.chdir(sub)
    ddGsim, ddGsim_err, M = get_ddG_matrix_M()
    ddG = ddGexp - ddGsim

    eps = model.contact_epsilons

    u,s,v = np.linalg.svd(M)
    s_norm = s/max(s)
    cutoffs = np.array(list(0.5*(s_norm[:-1] + s_norm[1:])) + [0.0])

    Mnorm = np.linalg.norm(M)
    cond_num = np.zeros(len(cutoffs),float)
    Xps = [] 
    Xp_cpxs = [] 
    ratios_xp = []
    ratios_cpx = []
    ## 2. If num_singular_values_include.txt exists
    ##   2a. Compute solution epsilon_prime and save in Mut_#/mut directory
    print "  Solving for new parameters. Parameters will be saved in mut/"
    print "# Sing. %10s %10s %10s" % ("Cond.Num.","|xp|/|eps|","|xp_cpx|/|eps|")
    solution_string = "# Sing. %10s %10s %10s\n" % ("Cond.Num.","|xp|/|eps|","|xp_cpx|/|eps|")
    for i in range(len(cutoffs)):

        cutoff = cutoffs[i]
        Mpinv = np.linalg.pinv(M,rcond=cutoff)
        Mpinvnorm = np.linalg.norm(Mpinv)
        cond_num[i] = Mnorm*Mpinvnorm

        x_particular = np.dot(Mpinv,ddG)
        np.savetxt("mut/xp%d.dat" % (i+1),x_particular)
        delta_eps_xp = x_particular
        ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)
        Xps.append(delta_eps_xp)
        ratios_xp.append(ratio_xp)

        try: 
            LP_problem, solution, x_particular, N = apply_constraints_with_cplex(model,ddG,M,cutoff)
            delta_eps = x_particular + np.dot(N,solution)
            ratio_cpx = np.linalg.norm(delta_eps)/np.linalg.norm(eps)
            Xp_cpxs.append(delta_eps)
            np.savetxt("mut/xp_cplex%d.dat" % (i+1),x_particular)
        except cplex.CplexSolverError:
            Xp_cpxs.append(np.zeros(len(x_particular)))
            ratio_cpx = 0.
        iteration_string = "%5d %10.5e %10.5f %10.5f" % (i+1,cond_num[i],ratio_xp,ratio_cpx)
        print iteration_string
        ratios_cpx.append(ratio_cpx)
        solution_string += iteration_string + "\n"
    open("mut/solution.log","w").write(solution_string) 
    plot_solution_info(model,cond_num,ratios_xp,ratios_cpx,Xps,Xp_cpxs)
    # Choose
    return Xps
    print "Succes"
    raise SystemExit

    ##   2b. Save model object with new contact parameters

    if not os.path.exists("mut/num_singular_values_include.txt"):
        print "ERROR!"
        print "  Need file mut/num_singular_values_include.txt", " to proceed"
        print "  Exiting"
        raise SystemExit
    else:
        temp = open("mut/num_singular_values_include.txt").read().rstrip("\n")
        if temp.endswith("xp"):
            compute_xp = True
            num_singular_values = int(temp.split()[0])
            savebeadbeadxp = "mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = savebeadbeadxp
        else:
            compute_xp = False
            num_singular_values = int(temp)
            savebeadbeadxp = "mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = "mut/"+newbeadbead
        cutoff = cutoffs[num_singular_values-1]
        print "  Using ",num_singular_values," singular values. Cutoff of = ",cutoff

    append_log(model.subdir,"Starting: Calculating_MC2004") 

    if compute_xp == True:
        ## Load solution for requested number of singular values.
        print "  Using ONLY the particular solution x_p as delta_eps!"
        Mpinv = np.linalg.pinv(M,rcond=cutoff)
        x_particular = np.dot(Mpinv,ddG)
        np.savetxt("mut/x_p.dat",x_particular)

        delta_eps = x_particular
        ratio = np.linalg.norm(delta_eps)/np.linalg.norm(eps)

        delta_eps_xp = x_particular
        ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)
    else:
        print "  Applying CPLEX to the particular solution to get delta_eps!"
        LP_problem, solution, x_particular, N = apply_constraints_with_cplex(model,ddG,M,cutoff)

        print "    Solution found!"
        delta_eps = x_particular + np.dot(N,solution)
        ratio = np.linalg.norm(delta_eps)/np.linalg.norm(eps)

        delta_eps_xp = x_particular
        ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)

    np.savetxt("mut/delta_eps.dat",delta_eps)
    np.savetxt("mut/delta_eps_xp.dat",delta_eps_xp)
    print "    Norm of perturbation, |deps|/|eps| = ", ratio
    print "    Norm of just x_p perturbation, |x_p|/|eps| = ", ratio_xp

    print "  Saving new parameters as: ",savebeadbead
    save_new_parameters(sub,eps,delta_eps,delta_eps_xp,savebeadbead, savebeadbeadxp, T)

    model.contact_energies = savebeadbead

    append_log(model.subdir,"Finished: Calculating_MC2004") 

def apply_constraints_with_cplex(model,ddG,M,cutoff):
    """ Construct and solve a linear/quadratic programming problem for new parameters.

    Description:

        Use the IMB ILOG Cplex optimization library to search the nullspace 
    dimensions for a solution that satisfies the given constraints and mimimizes
    the given objective function

    """

    eps = model.contact_epsilons

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
    Mpinv = np.linalg.pinv(M,rcond=cutoff)
    x_particular = np.dot(Mpinv,ddG)

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## S[:M.shape[1],:M.shape[1]] = np.diag(s)
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
    u,s,v = np.linalg.svd(M)
    rank = len(s)
    #print "  Matrix M singular values saved as: mut/singular_values.dat (and as ..._norm.dat)"
    #print "  Normed singular value spectrum:"
    #print s/max(s)

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

    ## Populate cplex linear (quadratic) programming problem
    LP_problem = cplex.Cplex()
    LP_problem.set_results_stream(None)
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
    LP_problem.objective.set_quadratic(objective_quadratic_coefficients)

    ## Let cplex do work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution_lambda = LP_problem.solution.get_values()

    ## Print cplex summary
    #print "Cplex summary:"
    #print "status: ",status
    #print "solution:",solution_lambda
    
    return LP_problem, solution_lambda, x_particular, N

def plot_solution_info(model,cond_num,ratios_xp,ratios_cpx,Xps,Xp_cpxs):
    """ Plot solution condition number and mutual covariance."""

    ## Plot condition number versus number of included singular values
    ## The condition number quantifies the intrinsic instability in inverting
    ## (solving) a given linear system.
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(range(1,len(cond_num)+1),cond_num,'b',label="$cond(M)$")
    ax1.plot([1],[cond_num[0]],'r',label="$log_{10}(cond(M))$")
    ax1.legend(loc=2)
    ax1.set_xlim(1,len(cond_num)+1)
    ax1.set_xlabel("# sing values")
    ax1.set_ylabel("$cond(M)$")
    ax1.set_title(model.subdir+" Condition number $cond(M) = ||M||\\cdot||M^+||$")
    ax2.plot(range(1,len(cond_num)+1),np.log10(cond_num),'r',label="$log(cond(M))$")
    ax2.set_ylabel("$log_{10}(cond(M))$")
    plt.savefig("mut/condition_number.pdf")
    np.savetxt("mut/cond_num.dat",cond_num)

    plt.figure()
    plt.plot(range(1,len(cond_num)+1),ratios_xp,color='r',marker='o',label="$|x_p|$/$|\\epsilon|$")
    plt.plot(range(1,len(cond_num)+1),ratios_cpx,color='b',marker='s',label="$|x_{cplex}|$/$|\\epsilon|$")
    plt.legend(loc=2)
    plt.title("Perturbation size $|\\delta\\epsilon|$/$|\\epsilon|$ w/ & w/o cplex")
    plt.savefig("mut/")
    temp = np.zeros((len(ratios_xp),2),float)
    temp[:,0] = np.array(ratios_xp)
    temp[:,1] = np.array(ratios_cpx)
    np.savetxt("mut/ratios.dat",temp)

    ## Compute correlation matrix for all possible xp's.
    ## Compute correlation matrix for all possible xp's with cplex.
    Covar_xp = np.zeros((len(Xps),len(Xps)),float)
    Covar_cpx = np.zeros((len(Xps),len(Xps)),float)
    Covar_xp_cpx = np.zeros((len(Xps),len(Xps)),float)
    for i in range(len(Xps)):

        xp1 = Xps[i]
        cp1 = Xp_cpxs[i]
        for j in range(len(Xps)):
            xp2 = Xps[j]
            cp2 = Xp_cpxs[j]

            Covar_xp[i,j]     = np.dot(xp1,xp2)/(np.linalg.norm(xp1)*np.linalg.norm(xp2))
            Covar_cpx[i,j]    = np.dot(cp1,cp2)/(np.linalg.norm(cp1)*np.linalg.norm(cp2))
            Covar_xp_cpx[i,j] = np.dot(xp1,cp2)/(np.linalg.norm(xp1)*np.linalg.norm(cp2))

    plt.figure()
    plt.pcolor(Covar_xp)
    plt.title("Covariance of $x_p$ and $x_p$")
    plt.colorbar()
    plt.savefig("mut/Covar_xp.pdf")
    plt.figure()
    plt.pcolor(Covar_cpx)
    plt.title("Covariance of $x_{cplex}$ and $x_{cplex}$")
    plt.colorbar()
    plt.savefig("mut/Covar_cpx.pdf")
    plt.figure()
    plt.pcolor(Covar_xp_cpx)
    plt.title("Covariance of $x_p$ and $x_{cplex}$")
    plt.colorbar()
    plt.savefig("mut/Covar_xp_cpx.pdf")

    plt.close()

def calculate_matrix_ddG_eps_M(model,coord):
    ''' Calculates and saves and returns the matrix from equation (9) in 
        Matysiak Clementi 2004. '''

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)

    os.chdir(model.subdir+"/mutants")
    print "  Loading mutants"
    mutants = mut.get_core_mutations()
    Fij, Fij_pairs, Fij_conts = phi.get_mutant_fij(model,mutants)

    os.chdir(sub)

    bounds, states = phi.get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

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

def get_ddG_matrix_M():
    ## Throw error if more than one temperature been used (?)
    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 
    if not os.path.exists("mut"):
        os.mkdir("mut")

    ## Collect ddG & M from Mut_#/<temp>_*/{mut,phi} directories.
    ## Take avg. and std error of mean for error bars. Save in Mut_#/mut

    print "  Getting experimental ddG"
    files = ["ddGsim.dat","ddGsim_err.dat","M.dat"] 
    flag = np.array([ not os.path.exists("mut/"+file) for file in files ])
    if np.any(flag):
        print "  One of the following does not exist: M.dat, ddG.dat. Calculating."
        ddGlist = []
        Mlist = []
        for n in range(len(directories)):
            dir = directories[n]
            ddG_temp = np.loadtxt(dir+"/phi/Q_phi.dat",usecols=(4,5),dtype=float)
            ddG = np.concatenate((ddG_temp[:,0],ddG_temp[:,1]),axis=0)
            ddGlist.append(ddG)
            M = np.loadtxt(dir+"/mut/M.dat",dtype=float)
            Mlist.append(M)
        #print "  IMPORTANT: Look at the singular value spectrum and choose a number of singular values to use."
        #print "  IMPORTANT: Make sure Mut_#/mut/num_singular_values_include.txt  exists."
        ddGlist = np.array(ddGlist)
        ddGsim = sum(ddGlist)/float(n+1)
        ddGsim_err = np.std(ddGlist,axis=0)/np.sqrt(float(n+1))
        Mlist = np.array(Mlist)
        M = sum(Mlist)/float(n+1)
        Merr = np.std(Mlist,axis=0)/np.sqrt(float(n+1))
        np.savetxt("mut/ddGsim.dat",ddGsim)
        np.savetxt("mut/ddGsim_err.dat",ddGsim_err)
        np.savetxt("mut/M.dat",M)
        np.savetxt("mut/Merr.dat",Merr)
    else:
        print "  Loading M.dat, ddG.dat, eps.dat"
        ddGsim= np.loadtxt("mut/ddGsim.dat")
        ddGsim_err = np.loadtxt("mut/ddGsim_err.dat")
        M = np.loadtxt("mut/M.dat")

    return ddGsim, ddGsim_err, M


if __name__ == '__main__':

    def dummy_func(sub,string):
        pass 
    
    subdirs = ["r15"]
    Models = models.load_models(subdirs,dryrun=True)
    model = Models[0]

    Xps = calculate_MC2004_perturbation(model,dummy_func)
    #calculate_matrix_ddG_eps_M(model,"Q")
