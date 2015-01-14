""" Compute different sets of native contacts

Description:

    This program calculates the native contacts of a trajectory (Gromacs .xtc
trajectory) given a BeadBead.dat file (which contains the equilibrium distances
for each non-bonded interaction) and Native.pdb (the native/crystal structure;
which holds the pairwise distances of the crystal structure).
    There are many ways to categorize contacts. For example, we can split
contacts into Native, Not-Native, Non-Native and sub-divide Native into
Native-Helical (secondary contacts) and Native-Non-Helical (tertiary contacts).


Procedure:
    Note: Must calculate the reference matrix first. 
    Note: The number of processors requested "mpirun -n 10" by match the number
    of temperatures to calculate.
1. First calculate the reference matrix at a given temperature (choose low T).
    python contact_calculator.py ref --method prob --refT ###
2. Submit a PBS job that executes the following command:
    mpirun -n 10 python contact_calculator.py  --method prob --refT ### --Ti ### --dT 10 --Tf ###


DEPRECATED AUG 2014

"""

import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

import mdtraj as md

def main():
    """ One possible branches: Calculate Q """
    parser = argparse.ArgumentParser(description='Calculate the (Non)Native contact matrix')

    parser.add_argument('--calc', action='store_true', help='calculate Q')

    args = parser.parse_args()
    
    if args.calc == True:
        calculate_Q()
    else:
        pass


def calculate_Q():
    """ Calculate contacts and save number of native contacts,
        number of native helical (and non-helical) contacts, 
        number of native local (and non-local) contacits. Uses
        MDTraj. WORKS!
    """
    
    Qref = np.loadtxt("Qref_cryst.dat")
    N = len(Qref)
    ## Turn contact matrix into vector by concatenating each row.
    native = []
    native_helical = []
    native_local = []
    for i in range(N-4):
        native.extend(Qref[i,i+4:])
        temp = list(np.zeros(len(Qref[i,i+4:])))
        temp2 = list(np.zeros(len(Qref[i,i+4:])))
        
        if len(temp) >= 5:
            temp[4] = 1
            temp2[4] = 1
            if len(temp2) >= 6:
                temp2[5] = 1
                if len(temp2) >= 7:
                    temp2[6] = 1
        native_helical.extend(temp)
        native_local.extend(temp2)
    native = np.array(native)
    native_helical = np.array(native_helical)
    native_local = np.array(native_local)

    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt("BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    print "  Computing distances with mdtraj..."
    traj = md.load("traj.xtc",top="Native.pdb")
    distances = md.compute_contacts(traj,pairs)
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)
    
    Qall = contacts*native
    Q = sum(Qall.T)
    A = sum((contacts*(1-native)).T)

    Qres = np.zeros((traj.n_frames,N),int) 
    Qhres = np.zeros((traj.n_frames,N),int) 
    Qlocalres = np.zeros((traj.n_frames,N),int) 
    print "  Summing native contacts per residue..."
    for k in range(N-4):
        slice = np.zeros(Qall.shape[1],int)
        sliceh = np.zeros(Qall.shape[1],int)
        slicelocal = np.zeros(Qall.shape[1],int)
        accum = 0
        for n in range(k+1):
            slice[accum+k-n] = 1
            
            if n == k:
                sliceh[accum+k-n] = 1
                slicelocal[accum+k-n] = 1
            elif (n == k-1) or (n == k-2):
                slicelocal[accum+k-n] = 1
            else:
                pass

            accum += N-4-n

        accum -= N-4-n
        Qres[:,k+4] = sum(Qall[:,slice==1].T)
        Qhres[:,k+4] = sum(Qall[:,sliceh==1].T)
        Qlocalres[:,k+4] = sum(Qall[:,slicelocal==1].T)

    Qnhres = Qres - Qhres  
    Qnonlocalres = Qres - Qlocalres 

    print "  Summing native contacts..."
    Qh = sum(Qhres.T)
    Qnh = sum(Qnhres.T)
    Qlocal = sum(Qlocalres.T)
    Qnonlocal = sum(Qnonlocalres.T)

    print "  Saving..."
    np.savetxt("Q.dat",Q)
    np.savetxt("A.dat",A)
    np.savetxt("Qh.dat",Qh)
    np.savetxt("Qnh.dat",Qnh)
    np.savetxt("Qlocal.dat",Qlocal)
    np.savetxt("Qnonlocal.dat",Qnonlocal)
    np.savetxt("Qres.dat",Qres,delimiter=" ",fmt="%d")
    np.savetxt("Qhres.dat",Qhres,delimiter=" ",fmt="%d")
    np.savetxt("Qnhres.dat",Qnhres,delimiter=" ",fmt="%d")
    np.savetxt("Qlocalres.dat",Qlocalres,delimiter=" ",fmt="%d")
    np.savetxt("Qnonlocalres.dat",Qnonlocalres,delimiter=" ",fmt="%d")

    ## Saving old filenames for backwards compatibility.
    np.savetxt("Qprob.dat",Q)
    np.savetxt("Qhprob.dat",Qh)
    np.savetxt("Qnhprob.dat",Qnh)

def get_beadbead_info(path='.'):
    """ Extract the native crystal structure contacts, equilibrium 
        contact distance (sigij), and number of residues N."""
    pairs = []
    #pdb = np.loadtxt("Native.pdb",dtype=str)
    #coords = pdb[:,6:9].astype(float)
    print "********* YOU SHOULDN'T SEE THIS MESSAGE!!! *************"
    coords = get_pdb_coords(path+"/Native.pdb")
    N = len(coords)
    Sig = np.zeros((N,N),float)
    Native_cryst = np.zeros((N,N),int)
    contacts = 0
    for line in open(path+"/BeadBead.dat","r"):
        i, j = int(line[0:5])-1, int(line[5:10])-1
        resi, resj = line[10:18], line[18:26]
        interaction_num = int(line[26:31])
        sigij, epsij = float(line[31:47]), float(line[47:63])
        Sig[i,j] = sigij
        rij = np.linalg.norm(coords[i] - coords[j])
        if rij <= 1.25*10.*sigij:
            Native_cryst[i,j] = 1

    #np.savetxt("Qref_cryst.dat",Native_cryst,fmt="%1d",delimiter=" ")
    return Native_cryst, Sig, N

def get_pdb_coords(pdbname):
    """ Get coordinates from a pdb file."""
    coords = []
    for line in open(pdbname,"r"):
        if line[:3] in ['TER','END']:
            break
        else:
            if line[:4] == "ATOM":
                coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 

    return np.array(coords)


if __name__ == "__main__":
    main()
