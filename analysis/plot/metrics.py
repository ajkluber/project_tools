import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import sys
import os
import shutil
import argparse
import numpy as np
import mdtraj as md


'''                                                                                                                               
Author: Fernando Yrazu                                                                                                            
Created: May 2014                                                                                                            
                                                                                                                                  
Description:                                                                                                                      
    This program plots several metrics (see description in main function of program) and saves the corresponding .pdf files. 
    It also saves the .dat files for inspection. This program should be executed from the parent directory of the simulation.
    It will store the files in a new directory called metrics within this parent directory. The Q_flow options will store 
    its output to /metrics/by_contact. The R_Q option will store its output to /metrics/rq

Procedure:

        1. Execute the program with the following options:                                                                                
        --type = (Required) Please input those metrics that you want to calculate (the options must be present in the dictionary list)     
        --prot = (Required) Input the proteins that you want to evaluate
        --iter = (Optional) Input the iteration number you are working with. Defaults to 0.

Changelog:                                                                                                                        
May 2014 Created                                                                                                             
June 2014 Added Leff                                                                                                                                  
'''


def get_args(metrics_dict, proteins_list):
    # User input to define which metrics are to be calculated, from a list displayed in the main() function
    
    parser = argparse.ArgumentParser(description='Select') 
    parser.add_argument('--type', type=str, required=True, choices=metrics_dict.keys(), nargs='+', help='Select metrics')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, nargs='+', help='Select proteins')
    parser.add_argument('--iter', type=str, required=False, help='Select iteration number')
    args = parser.parse_args()

    return args

def contacts_per_frame(protein, current_dir, select_path):
    beadbead = np.loadtxt(select_path+"BeadBead.dat",dtype=str)
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int)
    pairs -= np.ones(pairs.shape,int)
    np.savetxt(current_dir+'/metrics/by_contact/'+protein+"_contact_pairs.dat",pairs)

    print "Opening trajectory"
    traj = md.load(select_path+"traj.xtc",top=select_path+"Native.pdb")
    distances = md.compute_contacts(traj,pairs)
    all_contacts = (distances[0][:] <= 1.2*sigij).astype(int)
    local_contacts = (distances[0][:] <= 1.2*sigij).astype(int)
    non_local_contacts = (distances[0][:] <= 1.2*sigij).astype(int)
 
    print 'Contacts calculated for ', protein

   # Residue separation to be considered local:
    local_cutoff = 6 
    # Assign local and non_local contacts to the respective matrices
    for i in range(len(distances[1])):
        if all_contacts[0][i]==1:
            res_distance = abs(distances[1][i][0] - distances[1][i][1])
            # MDtraj already weeds out those residues less than (i, i+4) units apart
            if res_distance<=local_cutoff:
                local_contacts[0][i]=1
                non_local_contacts[0][i]=0
            elif res_distance>local_cutoff:
                local_contacts[0][i]=0
                non_local_contacts[0][i]=1                
    
    # The first row of every matrix contains the valid contacts. In the corresponding following functions, the rest of the rows
    # are trimmed correspondingly

    return all_contacts, local_contacts, non_local_contacts

def Q_avg(select_path, coordinate):

    Q_avg_frame = np.loadtxt(select_path+coordinate+'.dat')
    number_of_native_contacts = max(Q_avg_frame)
    Q_avg_frame = Q_avg_frame / number_of_native_contacts

    return Q_avg_frame, number_of_native_contacts

def partial_averages(Q_avg_frame, Q_joint):
    # Sort according to value of the first column (the value of Q)                                                                         
    col = 0 # corresponding to the value of Q per frame                                                                                    
    Q_sorted = Q_joint[Q_joint[:,col].argsort(axis=0)]
    Q_histogram, Q_bin_edges = np.histogram(Q_avg_frame, bins=50)

    #Find the average value of Q for each bin                                                                                              
    Q_midpoint = 0.5*(Q_bin_edges[1:]+Q_bin_edges[:-1])
    #Matrix to store the average Q value for each contact (Qi), for each bin                                                               
    Qi_matrix = np.zeros((len(Q_midpoint),len(Q_joint.T)-1))

    j = 0
    k = 0
    #Store the frame upper and lower indices in the sorted matrix, which correspond to the bin boundaries in the histogram               
    temparray = np.zeros(len(Q_bin_edges), dtype=int)

    for i in range(len(Q_bin_edges)-1):
        i+=1
        while Q_sorted[k,0] < Q_bin_edges[i]:
            k+=1

        temparray[j] = (k-1)
        j+=1

    # Final element in the list                                                                                                           
    temparray[j] = len(Q_joint)-1
    
    return Qi_matrix, Q_bin_edges, Q_sorted, Q_midpoint, temparray


def R_Q(protein, current_dir, select_path, iteration_number):
    print 'Calculating R_Q for ', protein
    # Calculate the R(Q) measure as per Chavez,Onuchic,Clementi(2004). For this we only need the Q evaluation per residue (Q_per_residue) 
    # as well as the overall Q_per_frame
    coordinate = 'Q'
    Q_avg_frame, number_of_native_contacts = Q_avg(select_path, coordinate)
  
    # Now proceed with the Q evaluation per residue. The Qlocal and Qnonlocal are not used in this function
    Q_per_contact_raw, Qlocal_per_contact_raw, Qnonlocal_per_contact_raw = contacts_per_frame(protein, current_dir, select_path)
    Q_per_contact = Q_avg_frame

    for i in range(len(Q_per_contact_raw.T)):
        if Q_per_contact_raw[0,i]==1:
            Q_per_contact = np.column_stack((Q_per_contact, Q_per_contact_raw[:,i]))
            
    Q_joint = Q_per_contact
    
    Qi_matrix, Q_bin_edges, Q_sorted, Q_midpoint, temparray = partial_averages(Q_avg_frame, Q_joint)
 
   # Initialize the matrix to store the R values 
    R_Q_single = np.zeros(len(Q_midpoint))
    
    # Calculate the R value for every Q_midpoint value

    for i in range(len(Q_bin_edges)-1):

        Qi_matrix[i,:] = np.average(Q_sorted[temparray[i]:temparray[i+1],1:], axis=0)
        summation = np.sum(np.square(Qi_matrix[i]-Q_midpoint[i]))/number_of_native_contacts 
        Q_product =  1/(Q_midpoint[i]*(1-Q_midpoint[i]))

        R_Q_single[i] = Q_product * summation
      
    R_Q = np.column_stack((Q_midpoint,R_Q_single))

    #Save the R_Q matrix
    np.savetxt(current_dir+'/metrics/rq/'+protein+'_'+iteration_number+'_R_Q.dat', R_Q)
        
    return R_Q 

def Leff(protein, current_dir, select_path, iteration_number):
    # We have the contact maps for the entire trajectory already calculated
    beadbead = np.loadtxt(select_path+"BeadBead.dat",dtype=str)
    pairs = beadbead[:,:2].astype(int)
    pairs -= np.ones(pairs.shape,int)
    native = beadbead[:,4].astype(int)

    pairs = pairs[ native != 0 ]

    number_of_contacts = len(pairs)
    # Raw sequence length for each contact
    sequence_dist = np.abs(pairs[:,1]-pairs[:,0])
    #Lj matrix of internal loops within contact loops
    lj = np.zeros((number_of_contacts, number_of_contacts))
    
    # Determine the raw number of internal loops within loops
    for i in range(len(lj)):
        pa = pairs[i][0]
        pb = pairs[i][1]
        for j in range(len(lj)):
            pa_bis = pairs[j][0]
            pb_bis = pairs[j][1]
            if any([all([pa_bis>=pa, pb_bis<pb]), all([pa_bis>pa, pb_bis<=pb])]):
                lj[i,j]=1
    
    # Weed out the intersecting loops
    for i in range(len(lj)):
        for j in range(len(lj)):
            if lj[i,j]==1:
                pa_bis = pairs[j][0]
                pb_bis = pairs[j][1]
                for k in range(len(lj)):
                    if lj[i,k]==1:
                        pa_tris = pairs[k][0]
                        pb_tris = pairs[k][1]
                        # If loop lj[i,k] is internal to lj[i,j], then set the internal loop to zero and continue
                        if (all([pa_tris>=pa_bis, pa_tris<pb_bis, pb_tris<pb_bis]) or 
                            all([pa_tris>pa_bis, pa_tris<pb_tris, pb_tris<=pb_bis])):
                            lj[i,k]=0
    # Not sure what to do with interlacing loops (there are very many)
    # Until the definitive answer is found:
    for i in range(len(lj)):
        for j in range(len(lj)):
            if lj[i,j]==1:
                pa_bis = pairs[j][0]
                pb_bis = pairs[j][1]
                for k in range(len(lj)):
                    if lj[i,k]==1:
                        pa_tris = pairs[k][0]
                        pb_tris = pairs[k][1]
                        if any([all([pa_tris<=pa_bis, pb_tris>pa_bis, pb_tris<pb_bis]),
                                all([pa_tris<pa_bis, pb_tris> pa_bis, pb_tris<=pb_bis]), 
                                all([pa_tris>pa_bis, pa_tris<pb_bis, pb_tris>=pb_bis]),
                                all([pa_tris>=pa_bis, pa_tris<pb_bis, pb_tris>pb_bis])]):
                            #Not taking into account probabilities here
                            if sequence_dist[j]>=sequence_dist[k]:
                                lj[i,k]=0
    # np.savetxt(current_dir+'/metrics/pCO_Leff/'+protein+'_lj.dat', lj, fmt='%3.0i',)
    lj_dist = lj * sequence_dist
    # Apply sequence distance vector to obtain true lj
    Q_avg_frame, number_of_native_contacts = Q_avg(select_path, 'Q')
    Q_per_contact_raw, Qlocal_per_contact_raw, Qnonlocal_per_contact_raw = contacts_per_frame(protein, current_dir, select_path)
    Q_per_contact = Q_avg_frame

    for i in range(len(Q_per_contact_raw.T)):
        if Q_per_contact_raw[0,i]==1:
            Q_per_contact = np.column_stack((Q_per_contact, Q_per_contact_raw[:,i]))

    Q_joint = Q_per_contact
    Qi_matrix, Q_bin_edges, Q_sorted, Q_midpoint, temparray = partial_averages(Q_avg_frame, Q_joint)

    for i in range(len(Q_bin_edges)-1):
        Qi_matrix[i,:] = np.average(Q_sorted[temparray[i]:temparray[i+1],1:], axis=0)

    L_eff = np.zeros(len(Q_midpoint))
    # Partial matrix of Lj * <Qj>
    L_partial = np.zeros(len(pairs))
    L_partial_Q1 = np.zeros(len(pairs))

    for i in range(len(Q_midpoint)):
        for j in range(len(pairs)):
            L_partial[j] = sequence_dist[j]-np.sum(lj_dist[j,:]*Qi_matrix[i,:])
    
        L_eff[i] = np.sum(L_partial)/np.sum(sequence_dist)

    # Same as above, only for Q=1 for calculating the pCO
    for i in range(len(pairs)):
        L_partial_Q1[i] = sequence_dist[i]-np.sum(lj_dist[i,:]*Qi_matrix[-1,:])
        L_eff_Q1 = np.sum(L_partial_Q1)/np.sum(sequence_dist)
   
    L_eff = np.column_stack((Q_midpoint, L_eff))
    return L_eff, L_eff_Q1, Qi_matrix, Q_midpoint

def pCO(protein, current_dir, select_path, iteration_number):

    L_eff, L_eff_Q1, Qi_matrix, Q_midpoint = Leff(protein, current_dir, select_path)
    number_of_contacts = float(len(Qi_matrix.T))
    p_CO = np.zeros(len(Q_midpoint))
    for i in range(len(Q_midpoint)):
        p_CO[i] = (1/number_of_contacts) * (L_eff[i,1]/L_eff_Q1) * np.sum(Qi_matrix[i,:])
    p_CO = np.column_stack((Q_midpoint, p_CO)) 
    return p_CO
    
def coordinate_contacts_flow(protein, current_dir, select_path, coordinate, iteration_number):
    
    # Frames are always indexed by Q (total native contacts). number_of_native_contacts is not used in this function
    Q_avg_frame, number_of_native_contacts = Q_avg(select_path, 'Q')

    # Only the select matrix is used
    Q_per_contact_raw, Qlocal_per_contact_raw, Qnonlocal_per_contact_raw = contacts_per_frame(protein, current_dir, select_path)
    # Should always be indexed by Q
    Qc_per_contact = Q_avg_frame
    
    if coordinate=='Q':
        Qc_per_contact_raw = Q_per_contact_raw
    elif coordinate=='Qlocal':
        Qc_per_contact_raw = Qlocal_per_contact_raw
    elif coordinate=='Qnonlocal':
        Qc_per_contact_raw = Qnonlocal_per_contact_raw
   
    # Select only native contacts for each frame
    for i in range(len(Qc_per_contact_raw.T)):
        if (Qc_per_contact_raw[0,i]==1):
            Qc_per_contact = np.column_stack((Qc_per_contact, Qc_per_contact_raw[:,i]))

    Q_joint = Qc_per_contact

    Qi_matrix, Q_bin_edges, Q_sorted, Q_midpoint, temparray = partial_averages(Q_avg_frame, Q_joint)

    for i in range(len(Q_bin_edges)-1):
        if temparray[i]==temparray[i+1]:
            pass
        else:            
            Qi_matrix[i,:] = np.average(Q_sorted[temparray[i]:temparray[i+1],1:], axis=0)

    Qc_flow = np.column_stack((Q_midpoint,Qi_matrix))
    
    return Qc_flow

def local_contacts_flow(protein, current_dir, select_path, iteration_number):
    coordinate = 'Qlocal'
    Qlocal_flow = coordinate_contacts_flow(protein, current_dir, select_path, coordinate, iteration_number)
    #Save the Qlocal matrix
    np.savetxt(current_dir+'/metrics/by_contact/'+protein+'_'+iteration_number+'_Qlocal_flow.dat', Qlocal_flow)

    return Qlocal_flow

def nonlocal_contacts_flow(protein, current_dir, select_path, iteration_number):
    coordinate = 'Qnonlocal'
    Qnonlocal_flow = coordinate_contacts_flow(protein, current_dir, select_path, coordinate, iteration_number)
    #Save the Qnonlocal matrix     
    np.savetxt(current_dir+'/metrics/by_contact/'+protein+'_'+iteration_number+'_Qnonlocal_flow.dat', Qnonlocal_flow)
    return Qnonlocal_flow

def total_contacts_flow(protein, current_dir, select_path, iteration_number):
    coordinate = 'Q'
    Q_flow = coordinate_contacts_flow(protein, current_dir, select_path, coordinate, iteration_number)
    #Save the Q matrix 
    np.savetxt(current_dir+'/metrics/by_contact/'+protein+'_'+iteration_number+'_Q_flow.dat', Q_flow)

    return Q_flow

def select_contacts_flow(*args):
    pass

def Q_pmf(protein, current_dir, select_path, iteration_number):
    shutil.copy2(select_path+'pmfs/Q_pmf.pdf', current_dir+'/metrics/pmfs/'+protein+'_iter_'+iteration_number+'_Q_pmf.pdf')

def plot_metrics(x, data_matrix, metrics_names, proteins_choice, current_dir, iteration_number):
    print 'Plotting for ', x
    # get the correct name of the protein 
    coordinate = metrics_names[x]
    print 'Chosen coordinate is: ', coordinate
    
    color_code = ['r-', 'g-', 'b-', 'm-', 'y-']
    plot_legend_y = {'R_Q':'R_Q','pCO':'pCO','Leff':'Leff','local_contacts_flow':'Contact index',
                     'nonlocal_contacts_flow':'Contact index',
                     'total_contacts_flow':'Contact index' ,'select_contacts_flow':'Contact index'}
    plot_legend_x = {'R_Q':'Q','pCO':'Q','Leff':'Q','local_contacts_flow':'Q', 'nonlocal_contacts_flow':'Q',
                     'total_contacts_flow':'Q' ,'select_contacts_flow':'Q'}

    if coordinate == 'R_Q':
        for i in proteins_choice:
            index = proteins_choice.index(i)
            R_Q = data_matrix[index]
            plt.plot(R_Q[1:,0],R_Q[1:,1],color_code[index], label=proteins_choice[index])
        
        plt.ylim((0,1))
        plt.xlabel(plot_legend_x[coordinate])
        plt.ylabel(plot_legend_y[coordinate])
        plt.legend()
        if len(proteins_choice)==1:
            proteins_text = proteins_choice[0]
            proteins_file = proteins_choice[0]
        elif len(proteins_choice)==2:
            proteins_text= proteins_choice[0]+' and '+proteins_choice[1]
            proteins_file = proteins_choice[0] + proteins_choice[1]
        elif len(proteins_choice)==3:
            proteins_text= proteins_choice[0]+', '+proteins_choice[1]+' and '+proteins_choice[2]
            proteins_file = proteins_choice[0] + proteins_choice[1] + proteins_choice[2]

        plt.title(coordinate + " for "+ proteins_text )
        plt.savefig(current_dir +"/metrics/rq/"+x+ "_" + proteins_file +"_"+iteration_number+".pdf")
        plt.clf()

    elif coordinate == "pCO":
        for i in proteins_choice:
            index = proteins_choice.index(i)
            p_CO = data_matrix[index]
            plt.plot(p_CO[:,0],p_CO[:,1],color_code[index], label=proteins_choice[index])

        plt.xlabel(plot_legend_x[coordinate])
        plt.ylabel(plot_legend_y[coordinate])
        plt.legend()
        if len(proteins_choice)==1:
            proteins_text = proteins_choice[0]
            proteins_file = proteins_choice[0]
        elif len(proteins_choice)==2:
            proteins_text= proteins_choice[0]+' and '+proteins_choice[1]
            proteins_file = proteins_choice[0] + proteins_choice[1]
        elif len(proteins_choice)==3:
            proteins_text= proteins_choice[0]+', '+proteins_choice[1]+' and '+proteins_choice[2]
            proteins_file = proteins_choice[0] + proteins_choice[1] + proteins_choice[2]

        plt.title(coordinate + " for "+ proteins_text )
        plt.savefig(current_dir +"/metrics/pCO_Leff/"+x+ "_" + proteins_file +"_"+iteration_number+".pdf")
        plt.clf()

    elif coordinate == "Leff":

        for i in proteins_choice:
            index = proteins_choice.index(i)
            L_eff = data_matrix[index][0]
            plt.plot(L_eff[:,0],L_eff[:,1],color_code[index], label=proteins_choice[index])

        plt.xlabel(plot_legend_x[coordinate])
        plt.ylabel(plot_legend_y[coordinate])
        plt.legend()
        if len(proteins_choice)==1:
            proteins_text = proteins_choice[0]
            proteins_file = proteins_choice[0]
        elif len(proteins_choice)==2:
            proteins_text= proteins_choice[0]+' and '+proteins_choice[1]
            proteins_file = proteins_choice[0] + proteins_choice[1]
        elif len(proteins_choice)==3:
            proteins_text= proteins_choice[0]+', '+proteins_choice[1]+' and '+proteins_choice[2]
            proteins_file = proteins_choice[0] + proteins_choice[1] + proteins_choice[2]

        plt.title(coordinate + " for "+ proteins_text )
        plt.savefig(current_dir +"/metrics/pCO_Leff/"+x+ "_" + proteins_file +"_"+iteration_number+".pdf")
        plt.clf()

    elif any([coordinate == "local_contacts_flow", coordinate=="nonlocal_contacts_flow", coordinate=="total_contacts_flow"]):
        for i in proteins_choice:
            Qc_flow = data_matrix[proteins_choice.index(i)]
            minQ = min(Qc_flow[:,0])
            maxQ = max(Qc_flow[:,0])
            incQ = (float(maxQ) - float(minQ))/len(Qc_flow)

            # You can switch off the interpolation by switching to interpolation='none'
            plt.pcolor(Qc_flow[:][1:].T, cmap='jet')
            cbar = plt.colorbar()
            cbar.set_label("Fraction individual contacts formed $Q_{individual}$")
            plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
            plt.ylabel(plot_legend_y[coordinate])
            plt.title(x + " for "+ i)
            plt.savefig(current_dir+'/metrics/by_contact/'+x+'_'+i+'_'+iteration_number+'.pdf')
            plt.clf()
            pass
    elif coordinate == "select_contacts_flow":
            pass

    
def main():
    
    # The dictionary of possible metrics in this module
    metrics_dict = {'R_Q':R_Q,'pCO':pCO,'Leff':Leff,'Qlocal':local_contacts_flow,
                    'Qnonlocal':nonlocal_contacts_flow, 
                    'Q':total_contacts_flow, 'select_contacts_flow':select_contacts_flow, 'Q_pmf': Q_pmf}

    metrics_names = {'R_Q':'R_Q','pCO':'pCO','Leff':'Leff','Qlocal':'local_contacts_flow',
                    'Qnonlocal':'nonlocal_contacts_flow',
                    'Q':'total_contacts_flow', 'select_contacts_flow':'select_contacts_flow', 'Q_pmf':'Q_pmf'}

    proteins_list = ['r15', 'r16', 'r17', '1SHG', '1RIS', '1TEN', '1K85','1E0G','1E41', 'sh3']
    
    # User selects which metrics to display based on the available options
    metrics_choice = [x for x in get_args(metrics_dict, proteins_list).type]
    proteins_choice = [x for x in get_args(metrics_dict, proteins_list).prot]

    iteration_number = get_args(metrics_dict, proteins_list).iter
    if iteration_number == None:
        iteration_number = '0'

    current_dir = os.getcwd()
    
    # We will create a separate "metrics" directory to store all the data
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')
    os.chdir(current_dir+'/metrics/')
    if os.path.isdir(current_dir+'/metrics/by_contact')==False:
        os.mkdir('by_contact')
    if os.path.isdir(current_dir+'/metrics/rq')==False:
        os.mkdir('rq')
    if os.path.isdir(current_dir+'/metrics/pCO_Leff')==False:
        os.mkdir('pCO_Leff')
    if os.path.isdir(current_dir+'/metrics/pmfs')==False:
        os.mkdir('pmfs')
    os.chdir(current_dir)

    #Using a rather convoluted way... not sure how to plug in the arguments of functions
    for x in metrics_choice:
        data_matrix = []
        if x == 'Q_pmf':
            for y in proteins_choice:
                select_temp = open(current_dir+'/'+y+'/Mut_'+iteration_number+'/T_array_last.txt').readline().split('_')[0]
                select_path = current_dir+'/'+y+'/Mut_'+iteration_number+'/'+select_temp+'_1/'
                f = metrics_dict.get(x)
                f(y,current_dir,select_path,iteration_number)
                
        else:
            f = metrics_dict.get(x)
            #Once the proper function is selected from the dictionary, the results are stored in the data matrix
            for y in proteins_choice:

                # Find out the proper path for each protein according to the folding temperature 
                select_temp = open(current_dir+'/'+y+'/Mut_'+iteration_number+'/T_array_last.txt').readline().split('_')[0]
                select_path = current_dir+'/'+y+'/Mut_'+iteration_number+'/'+select_temp+'_1/'

                data_matrix.append(f(y, current_dir, select_path,iteration_number))
            # And then feed into the plotting function
            plot_metrics(x,data_matrix,metrics_names, proteins_choice, current_dir,iteration_number)
    

if __name__ == "__main__":
    main()
















