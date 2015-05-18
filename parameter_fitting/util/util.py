'''Module containing all the shared function from the subdirectories of parameter_fitting

'''
import numpy as np
import time

try:
    import mdtraj as md
except:
    pass

def get_rij_Vp(model):
    ''' Load trajectory, state indicators, and contact energy '''
    #assumes you are in the directory with traj.xtc and Native.pdb
    time1 = time.time()
    traj = md.load("traj.xtc",top="Native.pdb")     # Loading from file takes most time.
    time2 = time.time()
    print " Loading traj took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)

    # rij is a matrix, where first index represents trajectory step, and the second index represents the different pairs
    time1 = time.time()
    rij = md.compute_distances(traj,model.pairs-np.ones(model.pairs.shape),periodic=False)
    time2 = time.time()
    print " Calculating rij took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)
    
    # NOT GENERALIZED FOR FITTING SUBSET OF MODEL PARAMETERS
    time1 = time.time()
    Vp = np.zeros((traj.n_frames,model.n_model_param),float)
    for i in range(model.n_pairs):   
        param_idx = model.pairwise_param_assignment[i]
        Vp[:,param_idx] = Vp[:,param_idx] + model.pair_V[i](rij[:,i])
    time2 = time.time()
    print " Calculating Vp took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)
    
    return traj,rij,Vp

def get_traj_rij(model):
    ''' Load trajectory and pairwise distances '''
    #assumes you are in the directory with traj.xtc and Native.pdb
    time1 = time.time()
    traj = md.load("traj.xtc",top="Native.pdb")     # Loading from file takes most time.
    time2 = time.time()
    print " Loading traj took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)

    # rij is a matrix, where first index represents trajectory step, and the second index represents the different pairs
    time1 = time.time()
    rij = md.compute_distances(traj,model.pairs-np.ones(model.pairs.shape),periodic=False)
    time2 = time.time()
    print " Calculating rij took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)
    
    return traj,rij

def get_rij(model,trajfiles,native):
    ''' Load trajectory and pairwise distances '''
    time1 = time.time()
    traj = md.load(trajfiles,top=native)     # Loading from file takes most time.
    time2 = time.time()
    print " Loading traj took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)

    # rij is a matrix, where first index represents trajectory step, and the second index represents the different pairs
    time1 = time.time()
    rij = md.compute_distances(traj,model.pairs-np.ones(model.pairs.shape),periodic=False)
    time2 = time.time()
    print " Calculating rij took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)
    
    return rij

def get_Vp_for_state(model,rij,state,n_frames):
    ''' Get contact energy for subset of frames'''
    
    time1 = time.time()
    Vp_state = np.zeros((n_frames,model.n_fitting_params),float)
    for i in range(model.n_fitting_params):   
        param_idx = model.fitting_params[i]

        # Loop over interactions that use this parameter
        for j in range(len(model.model_param_interactions[param_idx])):
            pair_idx = model.model_param_interactions[param_idx][j]
            Vp_state[:,i] = Vp_state[:,i] + model.pair_V[pair_idx](rij[state,pair_idx])
    time2 = time.time()
    print " Calculating Vp took: %.2f sec = %.2f min" % (time2-time1,(time2-time1)/60.)
    
    return Vp_state

def get_state_bounds():
    ''' Bounds for each state. Bounds are bin edges along Q. '''
    import os
    if not os.path.exists("state_bounds.txt"):
        raise IOError("Create state_bounds.txt with the boundaries of each state")
    else:
        statefile = open("state_bounds.txt","r").readlines()
    
    state_bounds = []
    state_labels = []
    for line in statefile:
        info = line.split()
        state_bounds.append(float(info[1]))
        state_bounds.append(float(info[2]))
        state_labels.append(info[0])
    
    return state_bounds,state_labels

def get_state_indicators(x,bounds):
    state_indicator = np.zeros(len(x),int)
    # Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (x > bounds[state_num]).astype(int)*(x <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(x))
    # Boolean arrays that indicate which state each frame is in.
    # States are defined by their boundaries along coordinate x.
    U  = ((x > bounds[1]).astype(int)*(x < bounds[2]).astype(int)).astype(bool)
    TS = ((x > bounds[3]).astype(int)*(x < bounds[4]).astype(int)).astype(bool)
    N  = ((x > bounds[5]).astype(int)*(x < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))

    return U,TS,N,Uframes,TSframes,Nframes

def concatenate_state_indicators(coorddirs,bounds,coord="Q.dat"):
    x = np.concatenate([ np.loadtxt("%s/%s" % (y,coord)) for y in coorddirs])
    state_indicator = np.zeros(len(x),int)
    # Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (x > bounds[state_num]).astype(int)*(x <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(x))
    # Boolean arrays that indicate which state each frame is in.
    # States are defined by their boundaries along coordinate x.
    U  = ((x > bounds[1]).astype(int)*(x < bounds[2]).astype(int)).astype(bool)
    TS = ((x > bounds[3]).astype(int)*(x < bounds[4]).astype(int)).astype(bool)
    N  = ((x > bounds[5]).astype(int)*(x < bounds[6]).astype(int)).astype(bool)
    Nframes  = sum(N.astype(int))
    Uframes  = sum(U.astype(int))
    TSframes = sum(TS.astype(int))

    return U,TS,N,Uframes,TSframes,Nframes
