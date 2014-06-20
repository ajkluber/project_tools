import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os
import argparse
import mdtraj as md
import sys


'''                                                                                                                                               
Author: Slight adaptation by Fernando Yrazu of code done mainly by Alex J Kluber                                                                                                                            
                                                                                                                                                  
Created: May 2014                                                                                                                                 
                                                                                                                                                  
Description:                                                                                                                                      
                                                                                                                                                  
    This program creates a modified BeadBead.dat file that has an extra column with the values corresponding to the intensity
    of the contact energies per residue. Thus, VMD can be used to visualize these energies using a color gradient.
    The output will be directed to /metrics/energy_vmd

Procedure:                                                                                                                                        
        1. Execute the program with the following options:                                                                                        
                                                                                                                                                  
        --prot = (Required) Input the proteins that you want to evaluate                                                                          
        --iter = (Optional) Input the iteration number you are working with. Defaults to 0.                                                       
                                                                                                                                                  
Changelog:                                                                                                                                        
                                                                                                                                                  
May 2014 Created                                                                                                                                  
                                                                                                                                                  
'''

def get_args(proteins_list):

    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, nargs='+', help='Select proteins')
    parser.add_argument('--iter', type=str, required=False, help='Select iteration number')
    args = parser.parse_args()

    return args


def modify_beadbead(current_dir, protein, iteration_number):
    path = current_dir + '/' + protein + '/Tf_'+iteration_number+'/'
    temp = open(path+'T_array.txt', 'r').readline().split()[0]
    path += temp +'/'

    beadbead = np.loadtxt(path+"BeadBead.dat",dtype=str)
    pairs = beadbead[:,:2].astype(int)
    pairs -= np.ones(pairs.shape)
    epsij = beadbead[:,6].astype(float)
    native = beadbead[:,4].astype(int)

    pairs = pairs[ native != 0 ]
    epsij = epsij[ native != 0 ]

    #Count the number of residues
    Qref = np.loadtxt(current_dir+'/'+protein+"/Qref_cryst.dat")
    N = len(Qref)
    C = np.zeros((N,N),float)
    
    for k in range(pairs.shape[0]):
        C[pairs[k,0],pairs[k,1]] = epsij[k]

    C_per_res = sum(C)
    C_per_res /= float(max(C_per_res))
        
    pdbfile = open(current_dir+'/'+protein+"/clean.pdb","r").readlines()

    newpdb = ""
    for i in range(len(pdbfile)):
        if pdbfile[i].startswith("END"):
            break
        else:
            oldline = pdbfile[i][:-1]
            #resnum = int(oldline[22:26]) 
            resnum = int(oldline[22:26]) - 1
            newline = oldline + ("%5.2f" % C_per_res[resnum]) + "\n"
            newpdb += newline
            #print len(oldline)
            #print oldline

    newpdb += "END"
    open(current_dir+"/metrics/energy_vmd/"+protein+"_residues_by_epsilons.pdb","w").write(newpdb)
    
if __name__ == "__main__":

    #List of possible proteins to choose from                                                                                     
    proteins_list = ['r15', 'r16', 'r17', '1SHG', '1RIS', '1TEN', '1K85','1E0G','1E41']
    proteins_choice = [x for x in get_args(proteins_list).prot]
    #Mutation number (defaults to 0)                                                                                              
    iteration_number = get_args(proteins_list).iter

    if iteration_number == None:
        iteration_number = '0'
    # Determine if the metrics directory exists and create it if not                                                              
    current_dir = os.getcwd()

    # We will create a separate "metrics/energy_vmd" directory to store all the data                                              
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')

    os.chdir(current_dir+'/metrics')
    if os.path.isdir(current_dir+'/metrics/energy_vmd')==False:
        os.mkdir('energy_vmd')
    os.chdir(current_dir)

    for protein in proteins_choice:
        modify_beadbead(current_dir, protein, iteration_number)

