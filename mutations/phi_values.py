import numpy as np
import os

import mdtraj as md


'''
Alexander Kluber

Seed script to figure out how to calculate Phi-values.
'''


def calculate_phi_values():

    modelname = 'wt'
    mutation_data = open("mutations.txt","r").readlines()[1:]
    mut_indx = [ mutation_data[i].split()[0] for i in range(len(mutation_data)) ]
    wt_res = [ mutation_data[i].split()[1] for i in range(len(mutation_data)) ]
    mut_res = [ mutation_data[i].split()[2] for i in range(len(mutation_data)) ]

    for i in range(len(mut_indx)):

        fij_kname = "fij_"+wt_res[i]+mut_indx[i]+mut_res[i]+".dat"
        print fij_kname
        fij_k = np.loadtxt(fij_kname)
    pass


def crunch_Hk():
    pass

def get_state_bounds(path,coord):
    ''' Get bounds for each state for specified coordinate.'''
    statefile = open(path+coord+"_states.txt","r").readlines()
    states = []
    for line in statefile:
        states.append([line.split()[0],float(line.split()[1]),float(line.split()[2])])
    return states
        
def load_eps_delta_sig_traj(subdir):
    ''' Load in the info from the BeadBead.dat file. Sig_ij, eps_ij, delta_ij and
        index pairs. This information is constant for a trajectory.'''
    print "Loading BeadBead.dat"
    beadbead = np.loadtxt(subdir+"BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    ## Use mdtraj to compute the distances between pairs.
    print "Loading traj.xtc"
    traj = md.load(subdir+"traj.xtc",top=subdir+"Native.pdb")
    print "Computing distances..."
    traj_dist = md.compute_distances(traj,pairs)

    ## Load in a sample fij matrix.
    print "Loading in fij_F90A.dat"
    fij_temp = np.loadtxt("r15/mutants/fij_F90A.dat")
    fij = []
    for i in range(len(fij_temp)-4):
        fij.extend(fij_temp[i,i+4:])
    fij = np.array(fij)

    return sigij,epsij,deltaij,pairs,traj,traj_dist

if __name__ == '__main__':
    ## First task is to calculate the perturbations for each mutation for
    ## each frame in the trajectory.

    R = 0.00831415
    T = "131.17"
    path = "r15/Mut_0/"

    if os.path.exists(path+T+"_phi") == False:
        os.mkdir(path+T+"_phi")

    subdir = path+T+"_1/"
    sigij,epsij,deltaij,pairs,traj,traj_dist = load_eps_delta_sig_traj(subdir)

    ## Loop over frames to calculate dH.
    print "computing dH and Qij"
    dH_k = np.zeros(traj.n_frames,float)
    #Qij = np.zeros(traj.n_frames,float)
    for j in range(traj.n_frames):
        rij = traj_dist[j]
        qij = 5.*((sigij/rij)**12) - 6.*deltaij*((sigij/rij)**10)
        dH_temp = sum(fij*epsij*qij)
        dH_k[j]  = dH_temp
        #Qij[j] = qij

    if os.path.exists("temp_dH.dat") == False:
        print "saving dH_k ..."
        np.savetxt("temp_dH.dat",dH_k)
    if os.path.exists("temp_Qij.dat") == False:
        print "saving Qij ..."
        #np.savetxt("temp_Qij.dat",Qij)
    print "Done."

