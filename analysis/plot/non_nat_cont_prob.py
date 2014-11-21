"""Visualization of the non-native contact probabilities for the different states

Description:
   This program outputs a non-native contact probability matrix for a selected protein and iteration 
   The format chosen is .png for storage space saving purposes.

Usage:
1. Execute the program with the following parameters:
   --prot = (Required) Input the protein that you want to use
   --iteration = (Required) Input the iteration number you are working with.

Changelog:
November 2014 Created by FMY
"""
import numpy as np
import mdtraj as md
import subprocess as sb
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
from itertools import izip

from project_tools.parameter_fitting.ddG_MC2004 import compute_Jacobian as cj
from model_builder import check_inputs as ci


def get_args():

    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True,  help='Select protein')
    parser.add_argument('--iteration', type=str, required=False, help='Select iteration number')
    args = parser.parse_args()

    return args

def find_non_native_contacts(native_pairs, num_nat_conts, traj):

    #pairs = now we have to calculate all the (i, i+4) and greater possibilities
    all_distances, all_pairs = md.compute_contacts(traj,contacts='all',scheme='ca')

    #all_pairs includes the i, i+3 contacts which we do not use in MC2004 
    #Thus, I have to weed out the (i,i+3) contacts
    counter=0
    for pair in all_pairs:
        if abs(pair[1]-pair[0])==3:
            counter+=1

    non_native_pairs = np.zeros((len(all_pairs)-len(native_pairs)-counter,2))
    
    #Now fill an array with the non-native pairs that we are using
    counter_nonnat = 0
    
    for pair in all_pairs:
        native = False
        if abs(pair[1]-pair[0])==3:
                native=True
        else:
            for nat_pair in native_pairs:
                if (pair[0]==nat_pair[0] and pair[1]==nat_pair[1]):
                    native=True
              
        if native==False:
            non_native_pairs[counter_nonnat,0]=pair[0]
            non_native_pairs[counter_nonnat,1]=pair[1]
            counter_nonnat+=1

    #Calculate average equilibrium distance between Calphas for native contacts and use it as
    # a parameter to compute non-native contacts
    avg_eq_dist = find_avg_eq_dist(native_pairs)
       
    non_native_distances, pairs = md.compute_contacts(traj, contacts=non_native_pairs, scheme='ca')
 
    return non_native_distances, non_native_pairs, avg_eq_dist
    
def find_avg_eq_dist(native_pairs):
    
    # Save the first frame as native.xtc to calculate the avg distance between Calphas
    if (os.path.exists('native.xtc')==False):
        sb.call('trjconv -f traj.xtc -o native.xtc -b 0 -e 0'.split())

    native_traj = md.load('native.xtc', top='Native.pdb')
    native_distances, pairs = md.compute_contacts(native_traj, contacts=native_pairs, scheme='ca')

    avg_eq_dist = np.average(native_distances)
 
    return avg_eq_dist

def find_non_native_cont_matrix(non_native_distances,avg_eq_dist,  U, TS, N, Uframes, TSframes, Nframes):
    #A contact is made when the residues are closer than 120% of the average of the native pairs eq. distance
    difference = 1.2*avg_eq_dist
    #First a boolean, then change to float so that multiplications can be made
    non_native_cont_matrix = non_native_distances<difference
    non_native_cont_matrix = non_native_cont_matrix.astype(float)
    
    #Each matrix is 500,001 by (#pairs), so I split the trajectory part in 3 to avoid memory issues
    non_native_cont_prob_U_1 = np.multiply(non_native_cont_matrix[:200000,:].T, U[:200000].astype(float))
    non_native_cont_prob_U_1 = np.sum(non_native_cont_prob_U_1, axis=1)
    non_native_cont_prob_U_2 = np.multiply(non_native_cont_matrix[200000:400000,:].T, U[200000:400000].astype(float))
    non_native_cont_prob_U_2 = np.sum(non_native_cont_prob_U_2, axis=1)
    non_native_cont_prob_U_3 = np.multiply(non_native_cont_matrix[400000:,:].T, U[400000:].astype(float))
    non_native_cont_prob_U_3 = np.sum(non_native_cont_prob_U_3, axis=1)
    #Merge count vectors
    non_native_cont_prob_U = non_native_cont_prob_U_1+non_native_cont_prob_U_2+ non_native_cont_prob_U_3
    #Divide by #frames to transform count into a probability
    non_native_cont_prob_U/=float(Uframes)

    non_native_cont_prob_TS_1 = np.multiply(non_native_cont_matrix[:200000,:].T, TS[:200000].astype(float))
    non_native_cont_prob_TS_1 = np.sum(non_native_cont_prob_TS_1, axis=1)
    non_native_cont_prob_TS_2 = np.multiply(non_native_cont_matrix[200000:400000,:].T, TS[200000:400000].astype(float))
    non_native_cont_prob_TS_2 = np.sum(non_native_cont_prob_TS_2, axis=1)
    non_native_cont_prob_TS_3 = np.multiply(non_native_cont_matrix[400000:,:].T, TS[400000:].astype(float))
    non_native_cont_prob_TS_3 = np.sum(non_native_cont_prob_TS_3, axis=1)
    non_native_cont_prob_TS = non_native_cont_prob_TS_1+non_native_cont_prob_TS_2+non_native_cont_prob_TS_3
    non_native_cont_prob_TS/=float(TSframes)
    
    non_native_cont_prob_N_1 = np.multiply(non_native_cont_matrix[:200000,:].T, N[:200000].astype(float))
    non_native_cont_prob_N_1 = np.sum(non_native_cont_prob_N_1, axis=1)
    non_native_cont_prob_N_2 = np.multiply(non_native_cont_matrix[200000:400000,:].T, N[200000:400000].astype(float))
    non_native_cont_prob_N_2 = np.sum(non_native_cont_prob_N_2, axis=1)
    non_native_cont_prob_N_3 = np.multiply(non_native_cont_matrix[400000:,:].T, N[400000:].astype(float))
    non_native_cont_prob_N_3 = np.sum(non_native_cont_prob_N_3, axis=1)
    non_native_cont_prob_N = non_native_cont_prob_N_1+non_native_cont_prob_N_2+non_native_cont_prob_N_3
    non_native_cont_prob_N/=float(Nframes)

    return non_native_cont_prob_U, non_native_cont_prob_TS, non_native_cont_prob_N

def plot_contact_probability(non_native_cont_prob_U, non_native_cont_prob_TS, non_native_cont_prob_N, non_native_pairs):
    number_of_res = np.max(non_native_pairs)+1
    non_native_cont_matrix_U = np.zeros((number_of_res, number_of_res), dtype=float)
    non_native_cont_matrix_TS = np.zeros((number_of_res, number_of_res), dtype=float)
    non_native_cont_matrix_N = np.zeros((number_of_res, number_of_res), dtype=float)

    #Fill the residue-residue contact matrix with the contact probability vector info
    for (pair,prob) in izip(non_native_pairs,non_native_cont_prob_U):
        i = pair[0]
        j = pair[1]
        non_native_cont_matrix_U[i,j]=prob
        non_native_cont_matrix_U[j,i]=prob

    #Plot for each state
    plt.figure()
    plt.pcolor(non_native_cont_matrix_U)
    plt.title('Non-native contact probability in U state')
    plt.xlim((0,number_of_res))
    plt.ylim((0,number_of_res))
    plt.colorbar()
    plt.xlabel('Residue #')
    plt.ylabel('Residue #')
    plt.savefig('non_nat_prob_U.png')
    plt.clf()

    for (pair,prob) in izip(non_native_pairs,non_native_cont_prob_TS):
        i = pair[0]
        j = pair[1]
        non_native_cont_matrix_TS[i,j]=prob
        non_native_cont_matrix_TS[j,i]=prob

    plt.figure()
    plt.pcolor(non_native_cont_matrix_TS)
    plt.title('Non-native contact probability in TS state')
    plt.xlim((0,number_of_res))
    plt.ylim((0,number_of_res))
    plt.colorbar()
    plt.xlabel('Residue #')
    plt.ylabel('Residue #')
    plt.savefig('non_nat_prob_TS.png')
    plt.clf()

    for (pair,prob) in izip(non_native_pairs,non_native_cont_prob_N):
        i = pair[0]
        j = pair[1]
        non_native_cont_matrix_N[i,j]=prob
        non_native_cont_matrix_N[j,i]=prob

    plt.figure()
    plt.pcolor(non_native_cont_matrix_N)
    plt.title('Non-native contact probability in N state')
    plt.xlim((0,number_of_res))
    plt.ylim((0,number_of_res))
    plt.colorbar()
    plt.xlabel('Residue #')
    plt.ylabel('Residue #')
    plt.savefig('non_nat_prob_N.png')
    plt.clf()



def assign_state(bounds,traj):
 
    #Use Q as the coordinate to determine states
    Q = np.loadtxt("Q.dat")
    state_indicator = np.zeros(len(Q),int)

    ## Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (Q > bounds[state_num]).astype(int)*(Q <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(Q))

    U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
    TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
    N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))

    return U,TS,N,Uframes,TSframes,Nframes

if __name__ == '__main__':
    

    curr_dir = os.getcwd()

    protein_choice = get_args().prot
    iteration_number = get_args().iteration
    
    #Extract the native contact pairs from model.info
    model = ci.load_model(curr_dir+'/'+protein_choice)
    native_pairs = model.contacts
    num_nat_conts = model.n_contacts
    # In mdtraj pairs are 0-indexed, while in model.contacts they are 1-indexed
    native_pairs -= [1,1]

    #Open the directory corresponding to the protein and iteration of choice
    os.chdir(protein_choice+'/Mut_'+iteration_number)

    #Take first directory of the three equilibrium runs 
    temp_dir = open('T_array_last.txt','r').readline().split()[0]
    
    #Determine state bounds
    state_bounds, state_labels = cj.get_state_bounds()
    
    #Load equilibrium trajectory
    os.chdir(temp_dir)    
    traj = md.load('traj.xtc', top='Native.pdb')
    
    #Determine non-native pairs, and their distances matrix across the trajectory.
    #Also, determine a contact whenever the distance between non-native residues is less than 120% of the 
    # average equilibrium distance of the native contacts
    non_native_distances, non_native_pairs, avg_eq_dist= find_non_native_contacts(native_pairs, num_nat_conts,traj)
    
    #Determine the frames belonging to each of the states (two-state in this case)
    bounds = [0] + state_bounds + [num_nat_conts]
    U,TS, N, Uframes, TSframes, Nframes = assign_state(bounds,traj)
    
    #Find out the probability of forming each non-native pair along the trajectory
    non_native_cont_prob_U, non_native_cont_prob_TS, non_native_cont_prob_N = find_non_native_cont_matrix(non_native_distances, avg_eq_dist, U, TS, N, Uframes, TSframes, Nframes)

    #Plot a heatmap of the non-native pair contact probability 
    plot_contact_probability(non_native_cont_prob_U, non_native_cont_prob_TS, non_native_cont_prob_N , non_native_pairs) 
    os.chdir(curr_dir)
