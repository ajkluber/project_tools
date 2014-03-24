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
    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt(subdir+"/BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    ## Use mdtraj to compute the distances between pairs.
    print "  Loading traj.xtc with mdtraj..."
    traj = md.load(subdir+"/traj.xtc",top=subdir+"/Native.pdb")
    print "  Computing distances with mdtraj..."
    traj_dist = md.compute_distances(traj,pairs)

    return sigij,epsij,deltaij,pairs,traj,traj_dist

def calculate_dH_for_mutants(Model,System,append_log):
    ''' First task is to calculate the perturbations for each mutation for
        each frame in the trajectory.'''

    
    R = 0.00831415
    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    if not os.path.exists(sub+"/Tf_choice.txt"):
        print "ERROR!"
        print "  Please create ",sub+"/Tf_choice.txt with your choice to do mutations at."
        print "  Exiting"
        raise SystemExit
    else:
        Tf_choice = open(sub+"/Tf_choice.txt").read().split()[0]
        print "  Calculating dH for temp ",Tf_choice

    T = Tf_choice
    savedir = sub+"/"+T+"_agg"
    os.chdir(System.subdir)

    mutants = [ x.split()[1]+x.split()[0]+x.split()[2] for x in open("mutants/mutations.txt","r").readlines()[1:] ]

    sigij,epsij,deltaij,pairs,traj,traj_dist = load_eps_delta_sig_traj(savedir)

    if Model.interaction_types[0] == "LJ12-10":
        def Qij(r,sig,delta):
            return 5.*((sig/r)**12) - 6.*delta*((sig/r)**10)
    else:
        print "ERROR!"
        print "  Unrecognized interaction type ", Model.interaction_types[0]
        print "  Exiting."
        raise SystemExit

    for mut in mutants:
        if not os.path.exists(savedir+"/dH_"+mut+".dat"):
            ## Load fij matrix.
            print "    Loading fij_"+mut+".dat"
            fij_temp = np.loadtxt("mutants/fij_"+mut+".dat")
            fij = []
            for i in range(len(fij_temp)-4):
                fij.extend(fij_temp[i,i+4:])
            fij = np.array(fij)

            ## Loop over frames to calculate dH.
            print "    Computing dH for ", mut
            dH_k = np.zeros(traj.n_frames,float)
            for j in range(traj.n_frames):
                rij = traj_dist[j]
                #qij = 5.*((sigij/rij)**12) - 6.*deltaij*((sigij/rij)**10)
                qij = Qij(rij,sigij,deltaij)
                dH_temp = sum(fij*epsij*qij)
                dH_k[j]  = dH_temp
            np.savetxt(savedir+"/dH_"+mut+".dat",dH_k)

    os.chdir(cwd)

if __name__ == '__main__':
    pass
