"""Module containing all the shared function from the subdirectories of parameter_fitting

"""
import mdtraj as md
import numpy as np

def get_rij_Vij(model):
    """ Load trajectory, state indicators, and contact energy """
    ##assumes you are in the directory with traj.xtc and Native.pdb
    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    ## rij is a matrix, where first index represents trajectory step, and the second index represents the different pairs
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape),periodic=False)
    
    Vp = np.zeros((traj.n_frames,model.n_contacts),float)
    for i in range(model.n_contacts):       ## loop over number of parameters
        param_idx = model.pairwise_param_assignment[i]
        Vp[:,param_idx] = Vp[:,param_idx] + model.pairwise_potentials[i](rij[:,i])
    Vij = Vp
    
    return traj,rij,Vij
