import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
import os
import mdtraj as md
import sys
import argparse

'''                                                                                                   
Author: Fernando Yrazu. Adaptation of code done by Alex J Kluber
                                                                                                                                  
Created: May 2014                                                                                                                 

Description:                                                                                          
                            
    This program plots the histograms and the stability of the proteins and saves the corresponding .pdf files.      
    It should be executed from the parent directory of the simulation.     
    It will store the files in a new directory called /metrics/stability within this parent directory.                                       
                                                                                                                                  
Procedure:                                                                                                                        
        1. Execute the program with the following options:

        --prot = (Required) Input the proteins that you want to evaluate                                                       
        --type = (Required) Please input those metrics that you want to calculate
        --iter = (Optional) Input the iteration number you are working with. Defaults to 0.                                                                                     
Changelog:                                                                                                                
         
May 2014 Created                                                                                                                  

'''
def get_args(metrics_dict, proteins_list):
    # User input to define which metrics are to be calculated, from a list displayed in the main() function                       

    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, nargs='+', help='Select proteins')
    parser.add_argument('--type', type=str, required=True, choices=metrics_dict.keys(), nargs='+', help='Select metrics')
    parser.add_argument('--iter', type=str, required=False, help='Select iteration number')
    args = parser.parse_args()

    return args

def determine_epsij(proteins_choice, iteration_number, current_dir):
    
    epsij_list = []
    local_epsij_list = []
    nonlocal_epsij_list = []
    # Get BeadBead.dat file from proteins
    for i in range(len(proteins_choice)):
        
        path = current_dir+'/'+proteins_choice[i]+'/Tf_'+iteration_number+'/'
        # All temperatures have the same BeadBead.dat file so I pick the first one from the T_array list 
        temp = open(path+'T_array.txt', 'r').readline().split()[0]
        # Store the BeadBead file in a list
        beadbead = np.loadtxt(path+temp+'/BeadBead.dat', dtype=str)
        pairs = beadbead[:,:2].astype(int)
        pairs -= np.ones(pairs.shape)
        epsij = beadbead[:,6].astype(float)
        native = beadbead[:,4].astype(int) 
        
        #Extract only the native contacts
        pairs = pairs[ native != 0 ]
        epsij = epsij[ native != 0 ]
        
        local_epsij = []
        nonlocal_epsij = []

        #Determine the locality or nonlocality of the contacts
        for i in range(len(pairs)):
            if abs(pairs[i,0]-pairs[i,1])<7:
                local_epsij.append(epsij[i])
            else:
                nonlocal_epsij.append(epsij[i])

        # Store the lists for each protein
        epsij_list.append(epsij)
        local_epsij_list.append(local_epsij)
        nonlocal_epsij_list.append(nonlocal_epsij)
    
    return epsij_list, local_epsij_list, nonlocal_epsij_list

def plot_hist_stab(proteins_choice, iteration_number, current_dir):
    
    epsij_list, local_epsij_list, nonlocal_epsij_list = determine_epsij(proteins_choice, iteration_number, current_dir)

    partial_name=''
    for i in proteins_choice:
        partial_name += i

    stab_file = open(current_dir+'/metrics/stability/stability_'+partial_name+'_iter_'+iteration_number+'.dat', 'w')

    for i in range(len(proteins_choice)):
        number_of_contacts = len(epsij_list[i])
        stability = np.sum(epsij_list[i])
        rel_stability = stability /number_of_contacts
        local_stability = np.sum(local_epsij_list[i])
        rel_local_stability = local_stability/number_of_contacts
        nonlocal_stability = np.sum(nonlocal_epsij_list[i])
        rel_nonlocal_stability = nonlocal_stability/number_of_contacts
        r_cd = nonlocal_stability/stability
        Qref = np.loadtxt(current_dir+'/'+proteins_choice[i]+"/Qref_cryst.dat")
        number_of_residues = len(Qref)
        tf = 36.081061 * stability/number_of_residues + 56.218196
        
        stab_file.write('Protein: '+proteins_choice[i]+'\n')
        stab_file.write('Number of contacts: %5d \n' % number_of_contacts)
        stab_file.write('Stability: %6.1f \n' % stability)
        stab_file.write('Relative stability: %6.3f \n' % rel_stability)
        stab_file.write('Local stability: %6.1f \n' % local_stability)
        stab_file.write('Relative local stability: %6.3f \n' % rel_local_stability)
        stab_file.write('Nonlocal stability: %6.1f \n' % nonlocal_stability)
        stab_file.write('Relative nonlocal stability: %6.3f \n' % rel_nonlocal_stability)
        stab_file.write('Ratio of local to total stability (R_c/d): %6.3f \n' % r_cd)
        stab_file.write('Expected Tf : %6.2f \n\n' % tf)
    stab_file.close()
    
    # Now let's plot the histograms
    color_code = ['r', 'g', 'b', 'm', 'y']

    plot_legend_x = 'Value of epsilon_ij '

    plt.figure()

    for i in range(len(proteins_choice)):
        plt.subplot(len(proteins_choice),1,i+1)
        plt.hist(epsij_list[i], bins=30, color=color_code[i], normed=True, label=proteins_choice[i])
        plt.legend()
        plt.ylim((0.0, 1.5))
    plt.xlabel(plot_legend_x)
    plt.suptitle("Histogram of contact energies for " + partial_name)
    plt.savefig(current_dir+'/metrics/stability/contact_energies_hist_'+partial_name+'.pdf')
    plt.clf()

def plot_energy_density_map(protein, iteration_number ,current_dir, select_path):

    colors = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]
    new_map = mpl.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

    colors_low = [('white')] + [(cm.winter(i)) for i in xrange(1,256)]
    low_map = mpl.colors.LinearSegmentedColormap.from_list('low_map', colors_low, N=256)

    colors_high = [('white')] + [(cm.summer(i)) for i in xrange(1,256)]
    high_map = mpl.colors.LinearSegmentedColormap.from_list('high_map', colors_high, N=256)
  
    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt(select_path+"BeadBead.dat",dtype=str)
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int)
    pairs -= np.ones(pairs.shape,int)
    native = beadbead[:,4].astype(int)

    pairs = pairs[ native != 0 ]
    epsij = epsij[ native != 0 ]
    pairs_low = pairs[ epsij<1. ]
    epsij_low = epsij[ epsij<1. ]
    pairs_high = pairs[ epsij>=1. ]
    epsij_high = epsij[ epsij>=1. ]
    
    Qref = np.loadtxt(current_dir+'/'+protein+"/Qref_cryst.dat")
    C = np.zeros(Qref.shape,float)
    C_low = np.zeros(Qref.shape,float)
    C_high = np.zeros(Qref.shape,float)

    print len(pairs_low), len(pairs_high)

    for k in range(len(pairs)):
        C[pairs[k][0],pairs[k][1]] = epsij[k]

    # C_low                                                                                                                      
    for k in range(len(pairs_low)):
        C_low[pairs_low[k][0],pairs_low[k][1]] = epsij_low[k]

    # C_high
    for k in range(len(pairs_high)):
        C_high[pairs_high[k][0],pairs_high[k][1]] = epsij_high[k]

    print "  Plotting..."
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    ax = plt.subplot(1,1,1,aspect=1)
    plt.pcolor(C,cmap=new_map, vmin=0, vmax=3)
    plt.xlim(0,len(Qref))
    plt.ylim(0,len(Qref))               
    ax = plt.gca()
    cbar = plt.colorbar()
#    cbar.set_clim(0,4)
    cbar.set_label("Contact energy density map",fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title('Contact energy density map for '+protein+', iteration '+iteration_number)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(15)
    print "  Saving..."
    plt.savefig(current_dir+'/metrics/stability/'+protein+"_contact_energy_density_map.pdf")
    plt.clf()

    print "  Plotting..."
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    ax = plt.subplot(1,1,1,aspect=1)
    plt.pcolor(C_low,cmap=low_map)
    plt.xlim(0,len(Qref))
    plt.ylim(0,len(Qref))
    ax = plt.gca()
    cbar = plt.colorbar()
    #    cbar.set_clim(0,4)                                                                                                           
    cbar.set_label("Contact energy density map",fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title('Contact energy density map for '+protein+', iteration '+iteration_number+', eps<1')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(15)
    print "  Saving..."
    plt.savefig(current_dir+'/metrics/stability/'+protein+"_low_contact_energy_density_map.pdf")
    plt.clf()
    
    print "  Plotting..."
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    ax = plt.subplot(1,1,1,aspect=1)
    plt.pcolor(C_high,cmap=high_map, vmin=1, vmax=3)
    plt.xlim(0,len(Qref))
    plt.ylim(0,len(Qref))
    ax = plt.gca()
    cbar = plt.colorbar()
    cbar.set_label("Contact energy density map",fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    plt.title('Contact energy density map for '+protein+', iteration '+iteration_number+', eps>=1')
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(15)
    print "  Saving..."
    plt.savefig(current_dir+'/metrics/stability/'+protein+"_high_contact_energy_density_map.pdf")
    plt.clf()




def main():
    #List of possible proteins and functions to choose from
    proteins_list = ['r15', 'r16', 'r17']
    metrics_dict = {'hist':plot_hist_stab, 'dens':plot_energy_density_map}

    proteins_choice = [x for x in get_args(metrics_dict, proteins_list).prot]
    metrics_choice = [x for x in get_args(metrics_dict, proteins_list).type]

   #Mutation number (defaults to 0)
    iteration_number = get_args(metrics_dict, proteins_list).iter
    if iteration_number == None:
        iteration_number = '0'

    # Determine if the metrics directory exists and create it if not
    current_dir = os.getcwd()

    # We will create a separate "metrics" directory to store all the data                                                        
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')

    os.chdir(current_dir+'/metrics/')
    if os.path.isdir(current_dir+'/metrics/stability')==False:
        os.mkdir('stability')
    os.chdir(current_dir)

    # Draw plots
    for x in metrics_choice:
        f = metrics_dict.get(x)
        if x =='dens':
            for y in proteins_choice:
                select_temp = open(current_dir+'/'+y+'/Tf_'+iteration_number+'/T_array.txt').readline().split()[0]
                select_path = current_dir+'/'+y+'/Tf_'+iteration_number+'/'+select_temp+'/'
                f(y, iteration_number, current_dir, select_path)

        else:
            f(proteins_choice,iteration_number, current_dir)
 
if __name__=='__main__':
    main()
