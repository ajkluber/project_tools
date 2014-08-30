""" Parameter fitting using Newton's method

Description:

    This submodule performs parameter fitting to match simulation features to
experimental (or designed) data. This parameter fitting is an extension of the
Matysiak, Clementi 2004 (see reference (1)). It uses the multivariate
Newton-Raphson method to correct simulation parameters in order to fit some
inputted target features. See references (2),(3) for more on Netwon's method in
general.


Submodules:

  newton 
      Solves for new parameters by Newton's method: by inverting the 
    Jacobian.

  ddG_MC2004
      Calculates Jacobian and feature vector for Delta Delta G's to match
    experimental data from phi-value analysis.


References:

1. Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.

2. J.E. Dennis; Robert B. Schnabel. "Numerical Methods for Unconstrained
Optimization and Nonlinear Equations". SIAM. 1996.

3. Jorge Nocedal; Stephen Wright. "Numerical Optimization". Springer. 2000


"""

import newton
import ddG_MC2004
import FRET
import RMSF
import contact_Qi

global available_methods
available_methods = ["ddG_MC2004","FRET","RMSF","contact_Qi"]

global modules
modules = {"ddG_MC2004":ddG_MC2004,"FRET":FRET,"RMSF":RMSF,"contact_Qi":contact_Qi}

def get_state_bounds():
    """ Bounds for each state. Bounds are bin edges along Q. """
    if os.path.exists("state_bounds.txt"):
        statefile = open("state_bounds.txt","r").readlines()
    else:
        print "ERROR!"
        print "  Please create state_bounds.txt"
        print "  With the boundaries of each state along Q"
        print "  Exiting"
        raise SystemExit
    
    state_bounds = []
    state_labels = []
    for line in statefile:
        info = line.split()
        state_bounds.append(float(info[1]))
        state_bounds.append(float(info[2]))
        state_labels.append(info[0])
    
    return state_bounds,state_labels

def get_states_Vij(model,bounds,epsilons,delta,sigmas):
    """ Load trajectory, state indicators, and contact energy """

    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape))
    Q = np.loadtxt("Q.dat")

    state_indicator = np.zeros(len(Q),int)
    ## Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (Q > bounds[state_num]).astype(int)*(Q <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(Q))
    ## Boolean arrays that indicate which state each frame is in.
    ## States are defined by their boundaries along coordinate Q.
    U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
    TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
    N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))

    ## Only count values of potential energy function where interaction is
    ## attractive.
    x = sigmas/rij
    x[(x > 1.09)] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero. 
    Vij = epsilons*(5.*(x**12) - 6.*deltas*(x**10))     ## To Do: Generalize to other contact functions

    return traj,U,TS,N,Uframes,TSframes,Nframes,Vij

def prepare_newtons_method(model,method,append_log):
    """ Prepare the files to do newtons method """

    if method not in available_methods:
        print "ERROR! requested method not found %s" % method
        print " Choose from available_methods: ",available_methods
        print " Exiting."
        raise SystemExit

    submodule = modules[method]

    name = model.subdir
    iteration = model.Mut_iteration

    append_log(name,"Starting: Calculating_Jacobian")

    target_feature, target_feature_err = submodule.compute_Jacobian.get_target_feature(model)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = submodule.compute_Jacobian.calculate_average_Jacobian(model)

    if not os.path.exists("%s/Mut_%d/newton" % (name,iteration)):
        os.mkdir("%s/Mut_%d/newton" % (name,iteration))
    if not os.path.exists("%s/Mut_%d/%s" % (name,iteration,method)):
        os.mkdir("%s/Mut_%d/%s" % (name,iteration,method))

    print "  Saving feature vector and Jacobian in %s/Mut_%d/%s" % (name,iteration,method)
    np.savetxt("%s/Mut_%d/%s/target_feature.dat" % (name,iteration,method), target_feature)
    np.savetxt("%s/Mut_%d/%s/target_feature_err.dat" % (name,iteration,method), target_feature_err)
    np.savetxt("%s/Mut_%d/%s/sim_feature.dat" % (name,iteration,method), sim_feature_avg)
    np.savetxt("%s/Mut_%d/%s/sim_feature_err.dat" % (name,iteration,method), sim_feature_err)
    np.savetxt("%s/Mut_%d/%s/Jacobian.dat" % (name,iteration,method), Jacobian_avg)
    np.savetxt("%s/Mut_%d/%s/Jacobian_err.dat" % (name,iteration,method) ,Jacobian_err)

    print "  Saving feature vector and Jacobian in %s/Mut_%d/newton" % (name,iteration)
    np.savetxt("%s/Mut_%d/newton/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/Mut_%d/newton/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/Mut_%d/newton/sim_feature.dat" % (name,iteration), sim_feature_avg)
    np.savetxt("%s/Mut_%d/newton/sim_feature_err.dat" % (name,iteration), sim_feature_err)
    np.savetxt("%s/Mut_%d/newton/Jacobian.dat" % (name,iteration), Jacobian_avg)
    np.savetxt("%s/Mut_%d/newton/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)

    append_log(name,"Finished: Calculating_Jacobian")
