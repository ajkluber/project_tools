'''Module containing all the shared function from the subdirectories of parameter_fitting

'''
import mdtraj as md
import numpy as np
import os

def get_rij_Vp(model):
    ''' Load trajectory, state indicators, and contact energy '''
    ##assumes you are in the directory with traj.xtc and Native.pdb
    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    ## rij is a matrix, where first index represents trajectory step, and the second index represents the different pairs
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape),periodic=False)
    
    Vp = np.zeros((traj.n_frames,model.n_model_param),float)
    for i in range(model.n_contacts):       ## loop over number of parameters
        param_idx = model.pairwise_param_assignment[i]
        Vp[:,param_idx] = Vp[:,param_idx] + model.pairwise_potentials[i](rij[:,i])
    
    return traj,rij,Vp

def get_state_bounds():
    ''' Bounds for each state. Bounds are bin edges along Q. '''
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

def get_state_indicators(x,bounds):
    state_indicator = np.zeros(len(x),int)
    ## Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (x > bounds[state_num]).astype(int)*(x <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(x))
    ## Boolean arrays that indicate which state each frame is in.
    ## States are defined by their boundaries along coordinate x.
    U  = ((x > bounds[1]).astype(int)*(x < bounds[2]).astype(int)).astype(bool)
    TS = ((x > bounds[3]).astype(int)*(x < bounds[4]).astype(int)).astype(bool)
    N  = ((x > bounds[5]).astype(int)*(x < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))

    return U,TS,N,Uframes,TSframes,Nframes
