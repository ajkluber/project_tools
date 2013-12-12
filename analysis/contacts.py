#!/usr/bin/env python
import numpy as np
import argparse
from coord_util import mol_reader
#import dmc_model.beadbead as bb
#from mpi4py import MPI

'''
Author: Alexander Kluber
Created: August 2013 

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

Changelog:
August 2013 Created
- Created to use as a command line utility.
October 2013
- Added contacts_for_states to calculate contact maps for a portion of a 
trajectory given the frame indices.
December 2013
** Copied over from dmc_model to model_builder. Changing around some stuff.

'''

def main():
    ''' Two possible branches: 1. Calculate reference matrix, 2. Calculate Q '''
    parser = argparse.ArgumentParser(description='Calculate the (Non)Native contact matrix')

    parser.add_argument('--calcQ', action='store_true', help='calculate Q')

    args = parser.parse_args()
    
    if args.calcQ == True:
        calculate_Q()
    else:
        pass
        


def contacts_for_states(framestate,numstates,cutoff=1.25):
    ''' Calculate contacts given a list of frames to use. This function is to be
        used by other programs. This function is still a little experimental.
    '''
    Native_cryst, Sig, N = get_beadbead_info()

    statesprobij = [ np.zeros((N,N),float) for i in range(numstates) ]
    #accum = [0,0,0]
    accum = np.zeros(numstates)
    frameidx = 0
    numframesused = 0
    frames = mol_reader.open("traj.xtc") 

    ## Loops over frames in trajectory.
    for frame in frames.read():
        state = framestate[frameidx] 
        if state == -1:
            pass
        else:
            D = np.ones((N,N),float)
            X = np.reshape(np.array(frame.coordinates), (N,3))
            for i in range(4,N):
                diff = X[0:-i,:] - X[i:,:]
                d = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
                ## Slice to insert values on the diagonal   
                diag = (np.arange(0,N-i),np.arange(i,N))
                D[diag] = d
            Contact = (D <= cutoff*10*Sig).astype(int)
            statesprobij[state] += Contact
            accum[state] += 1

        frameidx += 1
        if (frameidx % 20000) == 0:
            print "Finished"
            break
            print "Frame # ", frameidx
            if frameidx == 60000:
                break

    ## PROBABLY CHANGE TO JUST RETURN THE CONTACT MAPS
    for m in range(len(statesprobij)):
        print "Accum ", m , " ", accum[m], "   Length ", len(statesprobij[m]), " Counted ", framestate.count(m)
        if accum[m] != 0:
            statesprobij[m] /= float(accum[m])
        np.savetxt("contacts%s.dat" % m,statesprobij[m])
    return statesprobij

def calculate_Q():
    ''' A routine to calculate Q, Q_helical, Q_non-helical and A (not-native
        contacts) probabalistically. Meant to be run from a PBS script that 
        is submitted in the directory containing all the files.
    '''
    
    Native_cryst, Sig, N = get_beadbead_info()
    Native = np.loadtxt("Qref_prob.dat")
    ## Define Q_local as the helical contacts (i,i+4) as well as (i,i+5) 
    ## and (i,i+6). 
    h_diag = (np.arange(0,N-4),np.arange(4,N))
    #h_diag2 = (np.arange(0,N-5),np.arange(5,N))
    #h_diag3 = (np.arange(0,N-6),np.arange(6,N))
    Native_h = np.zeros((N,N),float)
    Native_h[h_diag] = Native[h_diag]
    #Native_h[h_diag2] = Native[h_diag2]
    #Native_h[h_diag3] = Native[h_diag3]

    Q = []
    Qh = []
    A = []
    framenum = 0 
    frames = mol_reader.open("traj.xtc") 
    for frame in frames.read():
        D = np.ones((N,N),float)
        X = np.reshape(np.array(frame.coordinates), (N,3))

        ## This loops allows the construction of the pairwise distance 
        ## distance matrix without a second inner loop by computing
        ## the diagonals. Takes advantage of vectorized operations.
        for i in range(4,N):
            diff = X[0:-i,:] - X[i:,:]
            d = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
            ## Slice to insert values on the diagonal   
            diag = (np.arange(0,N-i),np.arange(i,N))
            D[diag] = d

        ## Contact matrix. Distances are within 125% of equilibrium 
        ## distance. Sum up contacts to get Q, Q_helical, Q_non-helical.
        Contact = (D <= 1.25*10*Sig).astype(int)
        Q.append(sum(sum(Native*Contact)))
        A.append(sum(sum((1 - Native)*Contact)))
        Qh.append(sum(sum(Native_h*Contact)))

        framenum += 1
        #if (framenum % 1000) == 0:
        #    print "Frame #", framenum 

    np.savetxt("Qprob.dat",np.array(Q))
    np.savetxt("Qhprob.dat",np.array(Qh))
    np.savetxt("Qnhprob.dat",np.array(Q)-np.array(Qh))
    np.savetxt("Aprob.dat",np.array(A))
    return Q, Qh

def get_beadbead_info():
    ''' Extract the native crystal structure contacts, equilibrium 
        contact distance (sigij), and number of residues N.'''
    pairs = []
    #pdb = np.loadtxt("Native.pdb",dtype=str)
    #coords = pdb[:,6:9].astype(float)
    coords = get_pdb_coords("Native.pdb")
    N = len(coords)
    Sig = np.zeros((N,N),float)
    Native_cryst = np.zeros((N,N),int)
    for line in open("BeadBead.dat","r"):
        i, j = int(line[0:5])-1, int(line[5:10])-1
        resi, resj = line[10:18], line[18:26]
        interaction_num = int(line[26:31])
        sigij, epsij = float(line[31:47]), float(line[47:63])
        Sig[i,j] = sigij
        rij = np.linalg.norm(coords[i] - coords[j])
        if rij <= 1.25*10.*sigij:
            Native_cryst[i,j] = 1

    np.savetxt("Qref_cryst.dat",Native_cryst,fmt="%3d",delimiter=" ")
    return Native_cryst, Sig, N

def get_pdb_coords(pdbname):
    ''' Get coordinates from a pdb file.'''
    coords = []
    for line in open(pdbname,"r"):
        if line[:3] in ['TER','END']:
            break
        else:
            if line[:4] == "ATOM":
                coords.append([float(line[31:39]),float(line[39:47]),float(line[47:55])]) 

    return np.array(coords)

def probabilistic_reference(cutoff=1.25):
    ''' Calculate reference matrix probabilistically. Considers native 
        contacts as those that are formed in >=90% of the frames in 
        the native state ensemble. This is a subset of the crystal 
        structure contacts.

        The native ensemble is defined as all frames within 2.5 Angstrom
        RMSD to native (Note: the folded basin no longer has its minimum 
        at rmsd=0). Varying this cutoff between 2-3 Angstrom doesn't make
        a big difference. 
    '''
    Native_cryst, Sig, N = get_beadbead_info()
    numframes = -1
    numframesused = 0
    probij = np.zeros((N,N),float)

    frames = mol_reader.open("traj.xtc") 
    time, rmsd = np.loadtxt("rmsd.xvg",unpack=True)
    ## Loop over frames in trajectory.
    for frame in frames.read():
        numframes += 1
        if rmsd[numframes] > 0.25:
            continue
        else:
            D = np.ones((N,N),float)
            X = np.reshape(np.array(frame.coordinates), (N,3))
            ## This loops allows the construction of the pairwise distance 
            ## distance matrix without a second inner loop by computing
            ## the diagonals. Takes advantage of vectorized operations.
            for i in range(4,N):
                diff = X[0:-i,:] - X[i:,:]
                d = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
                ## Slice to insert values on the diagonal   
                diag = (np.arange(0,N-i),np.arange(i,N))
                D[diag] = d

            ## Contact matrix. Distances are within cutoff of equilibrium 
            ## distance.
            Contact= (D <= cutoff*10*Sig).astype(int)
            probij += Contact
            if (numframesused >= 1000) == 0:
                break
            numframesused+= 1
    Native = (probij >= .9*float(numframesused)).astype(int)
    np.savetxt("Qref_prob.dat",Native,delimiter=" ",fmt="%1d")
    return Native

if __name__ == "__main__":
    main()
