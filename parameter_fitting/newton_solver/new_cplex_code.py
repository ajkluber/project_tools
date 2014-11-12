
import matplotlib.pyplot as plt
import numpy as np
import os

import mdtraj as md
import cplex


#from project_tools.parameter_fitting.ddG_MC2004 import phi_values as phi
from project_tools.parameter_fitting.ddG_MC2004 import mutatepdbs as mut
from model_builder.models.SmogCalpha import SmogCalpha as sca

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
#    weight = 2.
#    ddGtarget = (ddGexp + (weight-1)*ddGsim)/2.
#    dg = ddGtarget - ddGsim
    dg = ddGexp - ddGsim
    np.savetxt("mut/ddGexp.dat",ddGexp)
#    np.savetxt("mut/ddGtarget.dat",ddGtarget)
    np.savetxt("mut/dg.dat",dg)
    np.savetxt("mut/ddGexp_err.dat",ddGexp_err)

    eps = model.contact_epsilons

    if not os.path.exists("mut/num_singular_values_include.txt"):
        #sca.append_log("Calculating MC2004, see")
        u,s,v = np.linalg.svd(M)
        s_norm = s/max(s)
        print "Singular values"
        for i in range(len(s)):
            print i,s[i]            
        print "Rank of M"
        print np.linalg.matrix_rank(M, tol=None)
        cutoffs = np.array(list(0.5*(s_norm[:-1] + s_norm[1:])) + [0.0])

        Mnorm = np.linalg.norm(M)
        cond_num = np.zeros(len(cutoffs),float)
        Xps = [] 
        Xp_cpxs = [] 
        ratios_xp = []
        ratios_cpx = []
        print "  Solving for new parameters. Parameters will be saved in mut/"
        print "# Sing. %10s %10s %10s %5s %10s" % ("Cond.Num.","|xp|/|eps|","|xp_cpx|/|eps|", "Status", "Cutoff")
        solution_string = "# Sing. %10s %10s %10s %5s %10s\n" % ("Cond.Num.","|xp|/|eps|","|xp_cpx|/|eps|", "Status", "Cutoff")
        
        max_solvable_eig = 0
        for i in range(len(cutoffs)):

            ## Generical solution uses pseudo-inverse of M.
            cutoff = cutoffs[i]
            Mpinv = np.linalg.pinv(M,rcond=cutoff)
            Mpinvnorm = np.linalg.norm(Mpinv)
            cond_num[i] = Mnorm*Mpinvnorm

            x_particular = np.dot(Mpinv,dg)
            print "Min of xp = ", np.min(x_particular)
            np.savetxt("mut/xp%d.dat" % (i+1),x_particular)
            delta_eps_xp = x_particular
            ratio_xp = np.linalg.norm(delta_eps_xp)/np.linalg.norm(eps)
            Xps.append(delta_eps_xp)
            ratios_xp.append(ratio_xp)
            plot_input_eigenvectors(v,i)
            plot_output_eigenvectors(u,i)
            ## Apply cplex
            try: 
                LP_problem, solution, x_particular_cpx, N, status, sensitivity = apply_constraints_with_cplex(model,dg,M,cutoff)

                if status ==1:
                    max_solvable_eig = i

                delta_eps = x_particular_cpx + np.dot(N,solution[:-1])
                ratio_cpx = np.linalg.norm(delta_eps)/np.linalg.norm(eps)
                Xp_cpxs.append(delta_eps)
                np.savetxt("mut/xp_cplex%d.dat" % (i+1),delta_eps)
            except cplex.exceptions.CplexSolverError:
                Xp_cpxs.append(np.zeros(len(x_particular)))
                ratio_cpx = 0.
            iteration_string = "%5d %10.5e %10.5f %10.5f %5d %10.5f" % (i+1,cond_num[i],ratio_xp,ratio_cpx,status, cutoff)
            print iteration_string
            print solution[-1]
            ratios_cpx.append(ratio_cpx)
            solution_string += iteration_string + "\n"
        open("mut/solution.log","w").write(solution_string) 
        plot_solution_info(model,s,cond_num,ratios_xp,ratios_cpx,Xps,Xp_cpxs)
    else:
        temp = open("mut/num_singular_values_include.txt").read().rstrip("\n")
        num_singular_values = temp.split()[0]
        print "  Using delta eps with %s singular values." % num_singular_values
        if temp.endswith("xp"):
            savebeadbeadxp = sub+"/mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = savebeadbeadxp
            delta_eps = np.loadtxt(sub+"/mut/xp"+num_singular_values+".dat")
        else:
            savebeadbeadxp = sub+"/mut/"+newbeadbead.split(".dat")[0]+"_xp.dat"
            savebeadbead = sub+"/mut/"+newbeadbead
            delta_eps = np.loadtxt(sub+"/mut/xp_cplex"+num_singular_values+".dat")

        print "  Saving new beadbead as: ",savebeadbead
        ## Write new beadbead to file
        model.contact_epsilons += delta_eps
#        model.contact_epsilons[model.contact_epsilons < 0.1] = 0.1
        model.get_pairs_string()
        open(savebeadbead,"w").write(model.beadbead)
        model.contact_energies = savebeadbead
        #sca.append_log("Finished: Calculating_MC2004") 

    plot_cumulative_solution(u,s,v,max_solvable_eig)
    os.chdir(cwd)

def apply_constraints_with_cplex(model,dg,M,cutoff):
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
    x_particular = np.dot(Mpinv,dg)

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## for i in range(len(s)):
    ##    S[i,i] = s[i]
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
    #objective_coeff = list(2.*np.dot(x_particular.T,N))

    ## Objective 3: Maximize attractiveness of minimal attraction interaction
    objective_coeff = np.hstack((np.zeros(N.shape[1], float), 1.0))

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
    # eps_lower_bound = 0.1
    # eps_upper_bound = 10.0
    # right_hand_side = list(eps + x_particular - eps_lower_bound*np.ones(len(eps)))
    # right_hand_side_2 = list(eps + x_particular - eps_upper_bound*np.ones(len(eps)))
    # right_hand_side.extend(right_hand_side_2)

    # column_names = [ "x"+str(i) for i in range(N.shape[1]) ] 

    # row_names = [ "c"+str(i) for i in range(2*N.shape[0]) ]
    # senses = "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)

    # rows = []
    # for n in range(2):
      #  for i in range(len(N)):
       #     rows.append([ column_names, list(-N[i,:]) ])
   
    ## Constraints version 5:                                                                                                                
    ## Require that (-EGap) be minimal.    

    ## (x_particular + x_n)_l +(eps)_l > eps_lower_bound
    ## --> x_n > -x_particular - eps + eps_lower_bound
    ## --> Sum(Lambda_i * Y_i)_l > -x_particular_l - eps_l + eps_lower_bound
    ##(x_particular + x_n)_l +(eps)_l < eps_upper_bound                                                                          
    ## --> x_n < -x_particular - eps + eps_upper_bound                                                                           
    ## --> Sum(Lambda_i * Y_i)_l < -x_particular_l - eps_l + eps_upper_bound  
    ## (x_particular + x_n)_l -EGap > 0.
    ## --> Sum(Lambda_i*Y_i)_l -EGap > -x_particular_l

    eps_lower_bound = 0.1
    eps_upper_bound = 10.0
    right_hand_side = list(-eps -x_particular + eps_lower_bound)
    right_hand_side_2 = list(-eps -x_particular + eps_upper_bound)
    right_hand_side.extend(right_hand_side_2)
    right_hand_side_3 = list(-x_particular)
    right_hand_side.extend(right_hand_side_3)

    column_names = [ "x"+str(i) for i in range(N.shape[1]) ]
    column_names.append("EGap")

    row_names = [ "c"+str(i) for i in range(len(right_hand_side)) ]
    senses = "G"*len(right_hand_side_2) + "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)

    rows = []

    for i in range(len(N)):
        temp = list(N[i,:])
        temp.append(float(0))
        rows.append([ column_names, temp ])
    for i in range(len(N)):
        temp = list(N[i,:])
        temp.append(float(0))
        rows.append([ column_names, temp ])
    for i in range(len(N)):
        temp = list(N[i,:])
        # The variable to be minimized is Z = (-EGap)
        temp.append(float(-1))
        rows.append([ column_names, temp ])

    ## Set quadratic terms in objective.
#    objective_quadratic_coefficients = [ 1. for j in range(N.shape[1]) ]

    ## Set upper and lower bounds on the solution. Arbitrary. Hopefully these 
    ## don't matter. These are bounds on vector lambda and the (-EGap)
    ## Since (-EGap) should be negative, the upper boundary for this variable is set to 0.
    upper_bounds = list(100.*np.ones(N.shape[1]+1))
#    upper_bounds.append(float(100))
    lower_bounds = list(-100.*np.ones(N.shape[1]+1))
#    lower_bounds.append(float(0))
#    print len(objective_coeff), len(upper_bounds), len(lower_bounds), len(column_names)
#    raise SystemExit
    ## Populate cplex linear (quadratic) programming problem
    LP_problem = cplex.Cplex()
    LP_problem.set_log_stream(None)
    LP_problem.set_results_stream(None)
#    LP_problem.set_error_stream('mut/e_'+str(cutoff)+'.log')
    LP_problem.objective.set_sense(LP_problem.objective.sense.maximize)
    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
#    LP_problem.objective.set_quadratic(objective_quadratic_coefficients)
#    LP_problem.objective.set_linear(objective_linear_coefficients)

#    print LP_problem.variables.get_names()
#    print LP_problem.variables.get_lower_bounds()
    ## Let cplex do work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    solution_lambda = LP_problem.solution.get_values()
    if status == 1:
        sensitivity = LP_problem.solution.sensitivity.objective("EGap")
    else:
        sensitivity = "N/A"
#    try:
#        conflict = LP_problem.conflict.get()
#    except: 
#        conflict = "No conflicts"
    ## Print cplex summary
    #print "Cplex summary:"
    #print "status: ",status
    #print "solution:",solution_lambda
    
    return LP_problem, solution_lambda, x_particular, N, status, sensitivity

def plot_input_eigenvectors(v, num_eigenvalue):
    i = num_eigenvalue
    plt.figure()
    plt.bar(range(len(v[i,:])),v[i,:])
    plt.title("Eigenvector components for eigenvalue #"+ str(i))
    plt.xlabel("Epsilon #")
    plt.savefig("mut/input_eig_"+str(i)+".pdf")
    plt.close()

def plot_output_eigenvectors(u,num_eigenvalue):
    i = num_eigenvalue
    plt.figure()
    plt.bar(range(len(u[:,i])),u[:,i])
    plt.title("Eigenvector components for eigenvalue #"+ str(i))
    plt.xlabel("Equation #")
    plt.savefig("mut/output_eig_"+str(i)+".pdf")
    plt.close()

def plot_cumulative_solution(u,s,v,max_solvable_eig):
    w = np.zeros(v.shape[1])
    x = np.zeros(v.shape[1])
    y = np.zeros(u.shape[0])
    z = np.zeros(u.shape[0])
    for i in range(max_solvable_eig):
        w += v[i,:]*s[i]
        x += abs(v[i,:]*s[i])
        y += u[:,i]*s[i]
        z += abs(u[:,i]*s[i])

    max_solvable_eig+=1

    plt.figure()
#    plt.bar(range(len(w)), w, color='r', label="Relative")
    plt.bar(range(len(x)), x, color='b', label="Absolute")
    plt.title("Cumulative influence of input for "+str(max_solvable_eig)+" eigenvalues")
    plt.xlabel("Epsilon #")
    plt.xlim((0,len(x)))
#    plt.legend()
    plt.savefig("mut/cumul_input_abs.pdf")
    plt.close()
    
    plt.figure()
    plt.bar(range(len(w)), w, color='r', label="Relative")
#    plt.bar(range(len(x)), x, color='b', label="Absolute")
    plt.title("Cumulative influence of input for "+str(max_solvable_eig)+" eigenvalues")
    plt.xlabel("Epsilon #")
    plt.xlim((0,len(w)))
#    plt.legend()
    plt.savefig("mut/cumul_input_rel.pdf")
    plt.close()

    plt.figure()
#    plt.bar(range(len(y)), y, color='r', label="Relative")
    plt.bar(range(len(z)), z, color='b', label="Absolute")
    plt.title("Cumulative influence of output for "+str(max_solvable_eig)+" eigenvalues")
    plt.xlabel("Constraint #")
    plt.xlim((0,len(z)))
#    plt.legend()
    plt.savefig("mut/cumul_output_abs.pdf")
    plt.close()

    plt.figure()
    plt.bar(range(len(y)), y, color='r', label="Relative")
#    plt.bar(range(len(z)), z, color='b', label="Absolute")
    plt.title("Cumulative influence of output for "+str(max_solvable_eig)+" eigenvalues")
    plt.xlabel("Constraint #")
    plt.xlim((0,len(y)))
#    plt.legend()
    plt.savefig("mut/cumul_output_rel.pdf")
    plt.close()

def plot_solution_info(model,s,cond_num,ratios_xp,ratios_cpx,Xps,Xp_cpxs):
    """ Plot solution condition number and mutual covariance."""

    plt.figure()
    plt.plot(s/max(s),'ro')
    plt.title(model.subdir+" Singular value spectrum for M")
    plt.savefig("mut/spectrum.pdf")
    np.savetxt("mut/singular_vals.dat",s)
    np.savetxt("mut/singular_vals_norm.dat",s/max(s))

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
    plt.xlabel("# sing values")
    plt.ylabel("Size of perturbation")
    plt.savefig("mut/perturb_size.pdf")
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
            nrmxp1 = np.linalg.norm(xp1)
            nrmxp2 = np.linalg.norm(xp2)
            nrmcp1 = np.linalg.norm(cp1)
            nrmcp2 = np.linalg.norm(cp2)

            if (nrmxp1 == 0.) or (nrmxp2 == 0.):
                Covar_xp[i,j]     = 0. 
            else:
                Covar_xp[i,j]     = np.dot(xp1,xp2)/(nrmxp1*nrmxp2)

            if (nrmcp1 == 0.) or (nrmcp2 == 0.):
                Covar_cpx[i,j]    = 0. 
            else:
                Covar_cpx[i,j]    = np.dot(cp1,cp2)/(nrmcp1*nrmcp2)

            if (nrmxp1 == 0.) or (nrmcp2 == 0.):
                Covar_xp_cpx[i,j]    = 0. 
            else:
                Covar_xp_cpx[i,j] = np.dot(xp1,cp2)/(nrmxp1*nrmcp2)

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

def save_new_parameters(sub,eps,delta_eps,delta_eps_xp,savebeadbead,savebeadbeadxp, T):
    """ Save new parameters as a BeadBead file

    DEPRECATED 7-7-14
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
    ## T_array_last.txt should hold names of directories 145.00_# where #=1,2,3
    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 
    if not os.path.exists("mut"):
        os.mkdir("mut")

    ## Collect ddG & M from Mut_#/<temp>_*/{mut,phi} directories.
    ## Take avg. and std error of mean for error bars. Save in Mut_#/mut
    print "  Getting simulation ddG"
    files = ["ddGsim.dat","ddGsim_err.dat","Jacobian.dat"] 
    flag = np.array([ not os.path.exists("mut/"+file) for file in files ])
    if np.any(flag):
        print "  One of the following does not exist: Jacobian.dat, ddG.dat. Calculating."
        ddGlist = []
        Mlist = []
        for n in range(len(directories)):
            dir = directories[n]
            ddG_temp = np.loadtxt(dir+"/phi/Q_phi.dat",usecols=(4,5),dtype=float)
            ddG = np.concatenate((ddG_temp[:,0],ddG_temp[:,1]),axis=0)
            ddGlist.append(ddG)
            M = np.loadtxt(dir+"/mut/Jacobian.dat",dtype=float)
            Mlist.append(M)
        #print "  IMPORTANT: Look at the singular value spectrum and choose a number of singular values to use."
        #print "  IMPORTANT: Make sure Mut_#/mut/num_singular_values_include.txt  exists."
        ddGlist = np.array(ddGlist)
        ddGsim = sum(ddGlist)/float(n+1)
        ddGsim_err = np.std(ddGlist,axis=0)
        Mlist = np.array(Mlist)
        M = sum(Mlist)/float(n+1)
        Merr = np.std(Mlist,axis=0)
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
