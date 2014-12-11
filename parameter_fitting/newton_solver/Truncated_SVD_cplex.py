
def find_solutions(model,method):
    target_feature = np.loadtxt("target_feature.dat")
    target_feature_err = np.loadtxt("target_feature_err.dat")
    sim_feature = np.loadtxt("sim_feature.dat")
    sim_feature_err = np.loadtxt("sim_feature_err.dat")
    Jacobian = np.loadtxt("Jacobian.dat")
    Jacobian_err = np.loadtxt("Jacobian_err.dat")
    J = Jacobian

    norm_eps = np.linalg.norm(model.model_param_values)

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
    
    for i in range(len(cutoffs)):
        S = np.zeros(J.shape) 
        tau = cutoffs[i] 
        s_use = s[s > tau]
        S[np.diag_indices(len(s_use))] = 1./s_use
        J_pinv = np.dot(v.T,np.dot(S.T,u.T))
        x_soln = np.dot(J_pinv,df)  ## particular
        
        #x_soln = np.dot(v.T,np.dot(S.T,np.dot(u.T,df)))
        try:
            LP_problem, solution, x_soln_cpx, status, sensitivity = apply_constraints_with_cplex(model,x_soln,N)
            if status == 1:
                max_solvable_eig = i
                residual = np.dot(J,x_soln) - df
                nrm_soln.append(np.linalg.norm(x_soln))
                nrm_resd.append(np.linalg.norm(residual))
                solutions.append(x_soln)
                Taus.append(tau)
        except:
            pass
    save_and_plot.save_solution_data(solutions,Taus,nrm_soln,nrm_resd,norm_eps,condition_number,s)

def apply_constraints_with_cplex(model,x_soln,N,weight=1.):
    """ Construct and solve a linear/quadratic programming problem for new parameters.

    Description:

        Use the IBM ILOG Cplex optimization library to search the nullspace 
    dimensions for a solution that satisfies the given constraints and mimimizes
    the given objective function

    """
    # Find correct wording for the number of native contacts
    num_native_contacts = model.n_native_contacts
    eps = model.model_param_values
    eps_native = eps[0:num_native_contacts]
    
    if len(eps)>len(eps_native):
        eps_non_native = eps[num_native_contacts:]
    else:
        eps_non_native = False


## ADDITION needed: non-native epsilons

    ## The general solution is a sum of the particular solution and an
    ## arbitrary vector from the nullspace of M.
#    Mpinv = np.linalg.pinv(M,rcond=cutoff)
    x_particular = x_soln

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

    
    if eps_non_native == False:
    ## Objective 3: Maximize attractiveness of minimal attraction interaction
        linear_objective_coeff = np.hstack((np.zeros(N.shape[1], float), -1.0))

    else:
    ## Objective_4: Maximize attractiveness of minimal native interaction (linear objective)
    ##              Minimize strength of non native interactions (quadratic objective)
    ## w is relative weith of frustration part with respect to native attractiveness part. We can only use a single
    ## objective with a cplex LP program. Thus, both objectives are simultaneous and a weight factor is added to 
    ## assign each of them a relative importance
    
        linear_objective_coeff = np.hstack((np.zeros(len(eps_native),float),np.zeros(len(eps_non_native),float), -1.0))
        quadratic_objective_coeff = np.hstack((np.zeros(len(eps_native),float),weight*np.ones(len(eps_non_native),float),0.0))

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
    eps_lower_bound_non_native = -3.
    eps_upper_bound_non_native = 3.

    if eps_non_native is False: 
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

    else:
        ## If we are considering the non native contacts, the restriction matrix changes
        right_hand_side = list(-eps_native -x_particular[0:num_native_contacts] + eps_lower_bound)
        right_hand_side_2 = list(-eps_native -x_particular[0:num_native_contacts] + eps_upper_bound)
        right_hand_side.extend(right_hand_side_2)
        right_hand_side_3 = list(-x_particular[0:num_native_contacts])
        right_hand_side.extend(right_hand_side_3)
        right_hand_side_4 = list(-eps_non_native -x_particular[num_native_contacts:] + eps_lower_bound_non_native)
        right_hand_side.extend(right_hand_side_4)
        right_hand_side_5 = list(-eps_non_native -x_particular[num_native_contacts:] + eps_upper_bound_non_native)
        right_hand_side.extend(right_hand_side_5)

        column_names = [ "x"+str(i) for i in range(len(eps)) ]
        column_names.append("EGap")

        row_names = [ "c"+str(i) for i in range(len(right_hand_side)) ]
        senses = "G"*len(right_hand_side_2) + "L"*len(right_hand_side_2) + "G"*len(right_hand_side_2)+"G"*len(right_hand_side_4)+"L"*len(right_hand_side_4)

        rows = []
        zeros_native = list(np.zeros(len(eps_native)))
        zeros_non_native = list(np.zeros(len(eps_non_native)))

        for i in range(len(N)):
            temp = list(N[i,:num_native_contacts])
            temp.extend(zeros_non_native)
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(N)):
            temp = list(N[i,:num_native_contacts])
            temp.extend(zeros_non_native)
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(N)):
            temp = list(N[i,:num_native_contacts])
            temp.append(zeros_non_native)
        # The variable to be minimized is Z = (-EGap)                                                                            
            temp.append(float(-1))
            rows.append([ column_names, temp ])
        for i in range(len(N)):
            temp = []
            temp.extend(zeros_native)
            temp2 = list(N[i,num_native_contacts:])
            temp.extend(temp2)
            temp.append(float(0))
            rows.append([ column_names, temp ])
        for i in range(len(N)):
            temp = []
            temp.extend(zeros_native)
            temp2 = list(N[i,num_native_contacts:])
            temp.extend(temp2)
            temp.append(float(0))
            rows.append([ column_names, temp ])
    
    
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
#    LP_problem.objective.set_sense(LP_problem.objective.sense.maximize)
    LP_problem.objective.set_sense(LP_problem.objective.sense.minimize)
#    LP_problem.variables.add(obj=objective_coeff, ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.variables.add(ub=upper_bounds, lb=lower_bounds, names=column_names)
    LP_problem.linear_constraints.add(lin_expr=rows, senses=senses, rhs=right_hand_side, names=row_names)
    LP_problem.objective.set_linear(linear_objective_coeff)
    if eps_non_native is not False:
        LP_problem.objective.set_quadratic(quadratic_objective_coeff)

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
    
    return LP_problem, solution_lambda, x_particular, status, sensitivity
