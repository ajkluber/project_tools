import numpy as np
import save_and_plot

def find_solutions(model,fitopts,position=100,eps_lower_bound=-2.,eps_upper_bound=3.,eps_lower_bound_non_native=-2.,eps_upper_bound_non_native=2.,EGap_lower_bound=-1.,frustration_fraction=0.):
    
    target_feature = np.loadtxt("target_feature.dat")
    target_feature_err = np.loadtxt("target_feature_err.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    sim_feature_err = np.loadtxt("sim_feature_err.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    Jacobian_err = np.loadtxt("Jacobian_err.dat")
    J = Jacobian

    model.n_long_native_pairs = model.n_long_pairs #Added
    binned_nonnatives = False
#    if J.shape[1] > (model.n_long_native_pairs):
#        binned_nonnatives = True
#        n_bins = J.shape[1]-model.n_long_native_pairs
    n_bins = J.shape[1]-model.n_long_native_pairs
    norm_eps = np.linalg.norm(model.long_model_param_values)

    ## Normalize the target step. 
    df = target_feature - sim_feature

    u,s,v = np.linalg.svd(J)
    N = v.T[:,J.shape[0]:]
    temp  = list(0.5*(s[:-1] + s[1:])) + [0.0]
    temp.reverse()
    cutoffs = np.array(temp)

    nrm_soln = []
    nrm_resd = []
    condition_number = []
    solutions = []
    Taus = []

    if (fitopts['nonnative']==True)==False or (binned_nonnatives == True):
        print "Solving for native parameters"
        for i in range(len(cutoffs)):                                                                                        
            S = np.zeros(J.shape)
            total_eigenvalues = len(list(cutoffs))
            n_eigenvalues = total_eigenvalues - i
            print "Position = "+str(i)+", taking "+str(n_eigenvalues)+" eigenvalues out of a total "+str(total_eigenvalues)
            tau = cutoffs[i]
            print "Cutoff value = "+str(tau)
            s_use = s[s > tau]
            S[np.diag_indices(len(s_use))] = 1./s_use
            J_pinv = np.dot(v.T,np.dot(S.T,u.T))
            x_particular = np.dot(J_pinv,df)  ## particular                                                                      
                                                                                               
            LP_problem,cplex_lambdas, status, sensitivity = apply_constraints_with_cplex(model,x_particular,N,eps_lower_bound,eps_upper_bound,eps_lower_bound_non_native,eps_upper_bound_non_native,EGap_lower_bound,frustration_fraction, binned_nonnatives, n_bins)
            print status
            if status == 1:
                break
    else:
        print "Solving for native and non native parameters"
        S = np.zeros(J.shape)
        i = position
        total_eigenvalues = len(list(cutoffs))
        n_eigenvalues = total_eigenvalues - position
        print "Position = "+str(i)+", taking "+str(n_eigenvalues)+" eigenvalues out of a total "+str(total_eigenvalues)
        tau = cutoffs[i] 
        print "Cutoff value = "+str(tau)
        s_use = s[s > tau]
        S[np.diag_indices(len(s_use))] = 1./s_use
        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_particular = np.dot(J_pinv,df)  ## particular
        
        LP_problem,cplex_lambdas, status, sensitivity = apply_constraints_with_cplex(model,x_particular,N,eps_lower_bound,eps_upper_bound,eps_lower_bound_non_native,eps_upper_bound_non_native,EGap_lower_bound,frustration_fraction, binned_nonnatives, n_bins)
        
    if status == 1:
        max_solvable_eig = i
        if len(cplex_lambdas)>(N.shape[1]+1):
            cplex_solution = x_particular + np.dot(N,cplex_lambdas[:-2])
            print "EGap = "+str(cplex_lambdas[-2])
            print "MaxFrustr" +str(-cplex_lambdas[-1])
        else:
            cplex_solution = x_particular + np.dot(N,cplex_lambdas[:-1])
            print "EGap = "+str(cplex_lambdas[-1])

        residual = np.dot(J,cplex_solution) - df
        nrm_soln.append(np.linalg.norm(cplex_solution))
        nrm_resd.append(np.linalg.norm(residual))
#        solutions.append(cplex_solution)
        Taus.append(tau)
        print "status:{0}, a feasible solution has been found".format(status)

        if binned_nonnatives == True:
            x_particular_binned = x_particular[:model.n_long_native_pairs]
            cplex_solution_binned = cplex_solution[:model.n_long_native_pairs]
            x_particular_reference = x_particular[model.n_long_native_pairs:]
            cplex_solution_reference = cplex_solution[model.n_long_native_pairs:]
            pairs_order = np.loadtxt('../../pairs_order.dat')

            for j in range(model.n_long_native_pairs,len(model.long_pairs)):
                order_index = pairs_order[j]
                x_particular_binned = np.append(x_particular_binned, x_particular_reference[order_index])
                cplex_solution_binned = np.append(cplex_solution_binned, cplex_solution_reference[order_index])
            
            x_particular = x_particular_binned
            cplex_solution = cplex_solution_binned

        solutions.append(cplex_solution)
        np.savetxt('x_particular.dat', x_particular)
        np.savetxt('lambda_vector.dat', cplex_lambdas)
        save_and_plot.save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s)                 
    else:
        print "status:{0} ".format(status)
        if status ==3:
            print "no feasible solution found"

#    np.savetxt('x_particular.dat', x_particular)
#    np.savetxt('lambda_vector.dat', cplex_lambdas)
#    save_and_plot.save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s)    
    parameters_log = open('cplex_parameters.dat','w')
    r = "Cplex simulation parameters\n\n"
    r+= "position = "+str(i)+"\n"
    r+= "eps_lower_bound = "+str(eps_lower_bound)+"\n"
    r+="eps_upper_bound = "+str(eps_upper_bound)+"\n"
    r+="eps_lower_bound_non_native = "+str(eps_lower_bound_non_native)+"\n"
    r+="eps_upper_bound_non_native = "+str(eps_upper_bound_non_native)+"\n"
    r+="EGap_lower_bound = "+str(EGap_lower_bound)+"\n"
    r+="frustration_fraction = "+str(frustration_fraction)+"\n"
    if status==1:
        r+= "Feasible solution found\n\n"
    else:
        r+= "No feasible solution found\n\n"
    parameters_log.write(r)
    open('Lambda_index.txt', 'w').write('0')
    open('desired_ratio', 'w').write('0.5')
    
def apply_constraints_with_cplex(model,x_particular,N,eps_lower_bound,eps_upper_bound,eps_lower_bound_non_native,eps_upper_bound_non_native,EGap_lower_bound,frustration_fraction, binned_nonnatives, n_bins):
    """ Construct and solve a linear/quadratic programming problem for new parameters.

    Description:

        Use the IBM ILOG Cplex optimization library to search the nullspace 
    dimensions for a solution that satisfies the given constraints and mimimizes
    the given objective function

    """
    import cplex

    # Find correct wording for the number of native contacts
    num_native_pairs = model.n_long_native_pairs
    eps = model.long_model_param_values
    eps_native = eps[0:num_native_pairs]
    
    if len(eps) == len(eps_native):
        eps_non_native = []
    elif binned_nonnatives == True:
        pairs_order = np.loadtxt('../../pairs_order.dat')
        eps_non_native = []
        for i in range(n_bins):
            for j in range(len(pairs_order[model.n_long_native_pairs:])):
                if pairs_order[j] == i:
                    eps_non_native.append(eps[j])
                    break
        eps_non_native = np.array(eps_non_native)
    else:
        eps_non_native = eps[num_native_pairs:]

#    print len(eps_non_native)
## ADDITION needed: non-native epsilons

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
#    Mpinv = np.linalg.pinv(M,rcond=cutoff)

    ## Singular value decomposition. As a test you can recover M by,
    ## S = np.zeros(M.shape)
    ## for i in range(len(s)):
    ##    S[i,i] = s[i]
    ## np.allclose(M,np.dot(u,np.dot(S,v))) --> should be True
#    u,s,v = np.linalg.svd(M)
#    rank = len(s)
    #print "  Matrix M singular values saved as: mut/singular_values.dat (and as ..._norm.dat)"
    #print "  Normed singular value spectrum:"
    #print s/max(s)

    ## Nullspace basis vectors are the last n-r columns of the matrix v.T. As a check
    ## all entries of the matrix np.dot(M,N) should be extremely small ~0. Because 
    ## they've been sent to the nullspace.
#    N = v.T[:,M.shape[0]:]

    ############# OBJECTIVE FUNCTION ###############
    ## Objective coefficients are sum of nullspace vectors. This comes from
    ## requiring the same average contact strength.

    ## Objective 1: Linear objective that seeks to minimize the change in 
    ## average contact strength. Found to give crazy results.
    #objective_coeff = list(sum(N))

    ## Objective 2: Quadratic objective that seeks to minimize the size of the
    ## perturbation.
    #objective_coeff = list(2.*np.dot(x_particular.T,N))

    
    if eps_non_native == []:
    ## Objective 3: Maximize attractiveness of minimal attraction interaction
        linear_objective_coeff = np.hstack((np.zeros(N.shape[1], float), 1.0))

    else:
    ## Objective_4: Maximize attractiveness of minimal native interaction (linear objective)
    ##              Minimize strength of non native interactions (quadratic objective)
    ## w is relative weith of frustration part with respect to native attractiveness part. We can only use a single
    ## objective with a cplex LP program. Thus, both objectives are simultaneous and a weight factor is added to 
    ## assign each of them a relative importance

#        linear_objective_coeff = np.hstack((np.zeros(N.shape[1],float),(np.ones(2))))
        linear_objective_coeff = np.hstack(( np.zeros(N.shape[1],float) , np.hstack((1-frustration_fraction,frustration_fraction)) ))
#        linear_objective_coeff = np.hstack((linear_objective_coeff, np.ones(len(eps_non_native),float)))
#        quadratic_objective_coeff = np.hstack((np.zeros(N.shape[1],float),-1.0))
#        quadratic_objective_coeff = np.hstack((quadratic_objective_coeff,np.ones(len(eps_non_native),float)))

    ############# CONSTRAINTS ###############f
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
    ## Require that (EGap) be maximal.    
    ## Require that (-Efrust) = NegEFrust be maximal. NegEfrust ranges from [-0.5,0]
                                           
    ## (x_particular + x_n)_l +(eps)_l > eps_lower_bound
    ## --> x_n > -x_particular - eps + eps_lower_bound
    ## --> Sum(Lambda_i * Y_i)_l > -x_particular_l - eps_l + eps_lower_bound
    ##(x_particular + x_n)_l +(eps)_l < eps_upper_bound                                                                          
    ## --> x_n < -x_particular - eps + eps_upper_bound                                                                           
    ## --> Sum(Lambda_i * Y_i)_l < -x_particular_l - eps_l + eps_upper_bound  
    ## (x_p + x_n)_l -EGap > 0.
    ## --> Sum(Lambda_i*Y_i)_l -EGap > -x_particular_l
    ## For the 4th condition (minimization of frustration):
    ## If contact was originally attractive:
    ## (x_p + x_n)_l + eps_non_native_l < Efrust = -NegEfrust
    ## --> x_n_l + NegEfrust < -x_p_l  - eps_non_native_l
    ## If contact was originally repulsive:
    ## (x_p + x_n)_l + eps_non_native_l > -Efrust = NegEfrust
    ## --> x_n_l - NegEfrust > -x_p_l  - eps_non_native_l  
    ## Enforce lower limits for originally attractive contacts and upper limits for originally repulsive contacts


    if eps_non_native == []: 
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
        # The variable to be maximized is EGap
            temp.append(float(-1))
            rows.append([ column_names, temp ])
    
    else:
        
        zeros_native = list(np.zeros(len(eps_native)))
        zeros_non_native = list(np.zeros(len(eps_non_native)))
        ## If we are considering the non native contacts, the restriction matrix changes
        right_hand_side = list(-eps_native -x_particular[0:num_native_pairs] + eps_lower_bound)
        right_hand_side_2 = list(-eps_native -x_particular[0:num_native_pairs] + eps_upper_bound)
        right_hand_side.extend(right_hand_side_2)
        right_hand_side_3 = list(-x_particular[0:num_native_pairs])
        right_hand_side.extend(right_hand_side_3)
        right_hand_side_4 = list(-eps_non_native -x_particular[num_native_pairs:])
        right_hand_side.extend(right_hand_side_4)
        right_hand_side_5 = list(-eps_non_native -x_particular[num_native_pairs:])
#        right_hand_side.extend(right_hand_side_5)
#        right_hand_side_6 = list(-eps_non_native -x_particular[num_native_pairs:]) 
#        right_hand_side.extend(right_hand_side_6)

        column_names = [ "x"+str(i) for i in range(N.shape[1])]
        column_names.append("EGap")
        column_names.append("NegEfrust")
#        column_names_2 = [ "x"+str(i) for i in range(len(eps_non_native))]
#        column_names.extend(column_names_2)

#        row_names = [ "c"+str(i) for i in range(len(right_hand_side)) ]
#        senses = "G"*len(right_hand_side_2) + "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)+"G"*len(right_hand_side_4)+"L"*len(right_hand_side_4) + "E"*len(right_hand_side_6)
#        senses = "G"*len(right_hand_side_2) + "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)+"L"*len(right_hand_side_4)
        senses = "G"*len(right_hand_side_2) + "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)
        rows = []

        for i in range(len(eps_native)):
            temp = list(N[i,:])
            temp.append(float(0))
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(eps_native)):
            temp = list(N[i,:])
            temp.append(float(0))
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(eps_native)):
            temp = list(N[i,:])
        # The variable to be maximized is EGap                                  
            temp.append(float(-1))
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(eps_native),len(N)):
            temp = list(N[i,:])
            temp.append(float(0))
            if model.long_pairwise_type[i] == 2:
                temp.append(float(1))
                senses = senses + "L"
            elif model.long_pairwise_type[i] == 1 or model.long_pairwise_type[i] == 3:
                temp.append(float(-1))
                senses = senses + "G"
            else:
                pass
            rows.append([ column_names, temp ])
            
        temp2 = []
        for i in range(len(eps_native),len(N)):
            temp = list(N[i,:])
            temp.append(float(0))
            temp.append(float(0))
            if model.long_pairwise_type[i] == 2:
                senses = senses + "G"
                temp2.append(-eps_upper_bound)
            elif model.long_pairwise_type[i] == 1 or model.long_pairwise_type[i] == 3:
                senses = senses + "L"
                temp2.append(eps_upper_bound)
            else:
                pass
            rows.append([ column_names, temp ])
        right_hand_side_5 = np.array(right_hand_side_5) + np.array(temp2)
        right_hand_side.extend(list(right_hand_side_5))

        row_names = [ "c"+str(i) for i in range(len(right_hand_side)) ]
#        for i in range(len(eps_native),len(N)):
#            temp = list(N[i,:])
#            temp.append(float(0))
#            temp.extend(zeros_non_native)
#            rows.append([ column_names, temp ])
        # constraints for minimization of eps_nonnative
#        for i in range(len(eps_non_native)):
#            temp = list(N[i,:])
#            temp.append(float(0))
#            temp.extend(zeros_non_native)
#            temp[N.shape[1]+1+i]=-1.
#            rows.append([ column_names, temp ])
            
#    print (np.array(rows)).shape
    ## Set upper and lower bounds on the solution. Arbitrary. Hopefully these 
    ## don't matter. These are bounds on vector lambda and the (-EGap)
    ## Since (-EGap) should be negative, the upper boundary for this variable is set to 0.
    upper_bounds = list(10000.*np.ones(N.shape[1]))
    upper_bounds.append(float(eps_upper_bound))
    if eps_non_native != []:
        upper_bounds.append(float(eps_upper_bound_non_native))
#    upper_bounds_2 = list(eps_upper_bound_non_native*np.ones(len(eps_non_native)))
#    upper_bounds.extend(upper_bounds_2)

    lower_bounds = list(-10000.*np.ones(N.shape[1]))
    lower_bounds.append(float(EGap_lower_bound))
    if eps_non_native != []:
        lower_bounds.append(float(eps_lower_bound_non_native))
#    lower_bounds_2 = list(eps_lower_bound_non_native*np.zeros(len(eps_non_native)))
#    lower_bounds.extend(lower_bounds_2)
#    print len(objective_coeff), len(upper_bounds), len(lower_bounds), len(column_names)
#    raise SystemExit
    ## Populate cplex linear (quadratic) programming problem
    LP_problem = cplex.Cplex()
#    LP_problem.set_log_stream(None)
#    LP_problem.set_results_stream(None)
#    LP_problem.set_error_stream('mut/e_'+str(cutoff)+'.log')
    LP_problem.objective.set_sense(LP_problem.objective.sense.maximize)
#    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
    print 'variables'
    LP_problem.variables.add(obj=linear_objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
#    LP_problem.variables.add(ub=upper_bounds, lb=lower_bounds, names=column_names)
    print 'constraints'
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
#    print 'objective'
#    linear_objective_pairs = zip(column_names, linear_objective_coeff)
#    LP_problem.objective.set_linear(linear_objective_coeff)
#    if eps_non_native == []:
#    LP_problem.objective.set_linear(linear_objective_pairs)
#    else:
#        quadratic_objective_pairs = zip(column_names, list(quadratic_objective_coeff))
#        LP_problem.objective.set_quadratic(quadratic_objective_pairs)
#        LP_problem.objective.set_linear(linear_objective_coeff)
#        LP_problem.objective.set_quadratic(quadratic_objective_coeff)
#    print LP_problem.variables.get_names()
#    print LP_problem.variables.get_lower_bounds()

    print 'solve'
    ## Let cplex do work.
    LP_problem.solve()
    status = LP_problem.solution.get_status()
    cplex_lambdas = LP_problem.solution.get_values()
    if eps_non_native == []:
        print "EGap" + str(cplex_lambdas[-1])
    else:
        print "EGap" + str(cplex_lambdas[-2])
        print "MaxFrust" + str(-cplex_lambdas[-1])

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
    
    return LP_problem, cplex_lambdas, status, sensitivity
