#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

from coord_util import mol_reader
import mdtraj as md


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
    ''' One possible branches: Calculate Q '''
    parser = argparse.ArgumentParser(description='Calculate the (Non)Native contact matrix')

    parser.add_argument('--calc', action='store_true', help='calculate Q')

    args = parser.parse_args()
    
    if args.calc == True:
        calculate_Q()
    else:
        pass

def contacts_for_states(framestate,numstates,name,cutoff=1.25):
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
            #print "Finished"
            #break
            print "Frame # ", frameidx
            #if frameidx == 60000:
            #    break

    ## PROBABLY CHANGE TO JUST RETURN THE CONTACT MAPS
    for m in range(len(statesprobij)):
        #print "Accum ", m , " ", accum[m], "   Length ", len(statesprobij[m]), " Counted ", list(framestate).count(m)
        if accum[m] != 0:
            statesprobij[m] /= float(accum[m])
        np.savetxt(name+"_%s.dat" % m,statesprobij[m])
    return statesprobij

def equil_contacts_for_states(framestate,numstates,savedir,Temp,numtemps,cutoff=1.25):
    ''' Calculate contacts given a list of frames to use. This function is to be
        used by other programs. This function is still a little experimental.
    '''

    Native_cryst, Sig, N = get_beadbead_info(path=Temp+"_1")

    statesprobij = [ np.zeros((N,N),float) for i in range(numstates) ]
    #accum = [0,0,0]
    accum = np.zeros(numstates)
    frameidx = 0
    numframesused = 0
    cwd = os.getcwd()
    ## Loop over subdirectories with same temp.
    for tnum in range(1,numtemps+1):
        T = Temp+"_"+str(tnum)
        
        os.chdir(T)

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
                #print "Finished"
                #break
                print "Frame # ", frameidx
                #if frameidx == 60000:
                #    break
        os.chdir(cwd)

    ## PROBABLY CHANGE TO JUST RETURN THE CONTACT MAPS
    for m in range(len(statesprobij)):
        #print "Accum ", m , " ", accum[m], "   Length ", len(statesprobij[m]), " Counted ", list(framestate).count(m)
        if accum[m] != 0:
            statesprobij[m] /= float(accum[m])
        np.savetxt(savedir+"/map_%s.dat" % m,statesprobij[m])
    return statesprobij

def calculate_Q():
    ''' Calculate contacts and save number of native contacts,
        number of native helical (and non-helical) contacts, 
        number of native local (and non-local) contacits. Uses
        MDTraj.
    '''
    
    Qref = np.loadtxt("Qref_cryst.dat")
    N = len(Qref)
    native = []
    native_helical = []
    native_local = []
    print len(Qref)
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
            sliceh[accum] = 1
            slicelocal[accum] = 1
            accum += N-4-n

        accum -= N-4-n
        Qres[:,k+4] = sum(Qall[:,slice==1].T)
        Qhres[:,k+4] = Qall[:,accum]
        if k >= N-6:
            if k < (N-5):
                Qlocalres[:,k+4] += Qall[:,accum]
                Qlocalres[:,k+4] += Qall[:,accum+1]
            else:
                Qlocalres[:,k+4] += Qall[:,accum]
        else:
            Qlocalres[:,k+4] += Qall[:,accum]
            Qlocalres[:,k+4] += Qall[:,accum+1]
            Qlocalres[:,k+4] += Qall[:,accum+2]

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

def old_calculate_Q():
    '''   DEPRECATED AK 4-2-2014
    '''
    
    Native_cryst, Sig, N = get_beadbead_info()
    #Native = np.loadtxt("Qref_prob.dat")
    Native = np.loadtxt("Qref_cryst.dat")
    ## Define Q_local as the helical contacts (i,i+4) as well as (i,i+5) 
    ## [and (i,i+6)]. This is assuming that contacts to i+4 and i+5 
    ## stabilize helices.
    h_diag = (np.arange(0,N-4),np.arange(4,N))
    h_diag2 = (np.arange(0,N-5),np.arange(5,N))
    #h_diag3 = (np.arange(0,N-6),np.arange(6,N))
    Native_h = np.zeros((N,N),float)
    Native_hi5 = np.zeros((N,N),float)
    Native_h[h_diag] = Native[h_diag]
    #Native_h[h_diag2] = Native[h_diag2]
    #Native_h[h_diag3] = Native[h_diag3]
    Native_hi5[h_diag2] = Native[h_diag2]

    Q = []
    Qh = []
    Qres = []
    Qhres = []
    Qhi5res = []
    A = []
    framenum = 0 
    frames = mol_reader.open("traj.xtc") 
    for frame in frames.read():
        D = np.zeros((N,N),float)
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
        Contact1 = (D <= 1.25*10*Sig).astype(int)
        Contact2 = (D > 0.01).astype(int)
        Contact = Contact1*Contact2
        Q.append(sum(sum(Native*Contact)))
        Qres.append(sum(Native*Contact))
        Qhres.append(sum(Native_h*Contact))
        Qhi5res.append(sum(Native_hi5*Contact))
        A.append(sum(sum((1 - Native)*Contact)))
        Qh.append(sum(sum(Native_h*Contact)))

        framenum += 1
        #if (framenum % 1000) == 0:
        #    print "Frame #", framenum 

    ## Save all the data files!
    np.savetxt("Qprob.dat",np.array(Q))
    np.savetxt("Qres.dat",np.array(Qres),delimiter=" ",fmt="%d")
    np.savetxt("Qhres.dat",np.array(Qhres),delimiter=" ",fmt="%d")
    np.savetxt("Qhi5res.dat",np.array(Qhi5res),delimiter=" ",fmt="%d")
    np.savetxt("Qhprob.dat",np.array(Qh))
    np.savetxt("Qnhprob.dat",np.array(Q)-np.array(Qh))
    np.savetxt("Aprob.dat",np.array(A))
    return Q, Qh

def get_beadbead_info(path='.'):
    ''' Extract the native crystal structure contacts, equilibrium 
        contact distance (sigij), and number of residues N.'''
    pairs = []
    #pdb = np.loadtxt("Native.pdb",dtype=str)
    #coords = pdb[:,6:9].astype(float)
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
        
        **Currently not used models based on Go model. 1-23-14 AK**

        coord_util mol_reader DEPRECATED 4-2-2014 AK

    '''
    Native_cryst, Sig, N = get_beadbead_info()
    numframes = -1
    numframesused = 0
    probij = np.zeros((N,N),float)
    print "## DEBUGGING: Crystal contacts", sum(sum(Native_cryst))
    print "## DEBUGGING: N = ", N
    print "## DEBUGGING: Initial probij", sum(sum(probij))

    frames = mol_reader.open("traj.xtc") 
    time, rmsd = np.loadtxt("rmsd.xvg",unpack=True)
    ## Loop over frames in trajectory.
    for frame in frames.read():
        numframes += 1
        if rmsd[numframes] > 0.25:
            continue
        else:
            #D = np.ones((N,N),float)  ## Why is D = ones((N,N))?
            D = np.zeros((N,N),float)
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
            Contact1 = (D <= 1.25*10*Sig).astype(int)
            Contact2 = (D > 0.01).astype(int)
            Contact = Contact1*Contact2
            #Contact = (D <= cutoff*10.*Sig).astype(int)

            probij += Contact
            if numframesused >= 1000:
                break
            numframesused+= 1

    probij /= float(numframesused)
    Native = (probij >= .9).astype(int)
    #print Native
    print "## DEBUGGING: Final probij", sum(sum(probij))
    print "## DEBUGGING: Probabalistic contacts: ",sum(sum(Native))
    print "## DEBUGGING: Type native",type(Native)
    print "## DEBUGGING: Native shape ", Native.shape
    #np.savetxt("probij.dat",probij,delimiter=" ")
    np.savetxt("Qref_prob.dat",Native,delimiter=" ")

    plt.pcolor(Native_cryst + probij.T)
    plt.xlim(0,len(Native_cryst))
    plt.ylim(0,len(Native_cryst))
    plt.ylabel("Contact Probability")
    plt.xlabel("Crystal Contacts")
    plt.title("1000 frames within 2.5A of Crystal")
    plt.savefig("crystal_vs_prob.pdf")
    
    return Native,probij

if __name__ == "__main__":
    main()
