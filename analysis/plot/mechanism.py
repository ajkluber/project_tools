""" Different ways to visualize the folding mechanism

Description:

    There are many ways to visualize the folding mechanism. 
    This programs saves its output files to a directory called /metrics/thermodynamic

Usage :
1. Execute the program with the following options:                                                               
        --prot = (Required) Input the proteins that you want to evaluate                                                              
        --iter = (Optional) Input the iteration number you are working with. Defaults to 0. 

Changelog:
May 2014 Created by AJK
May 2014 Modified by FY: Added Qlocal, Qnonlocal to plots
July 2014 Updated by FY to reflect changes in package directory structure
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import argparse
import mdtraj as md

def get_args(proteins_list):
                                                                                                                             
    parser = argparse.ArgumentParser(description='Select')
    parser.add_argument('--prot', type=str, required=True, choices=proteins_list, nargs='+', help='Select proteins')
    parser.add_argument('--iter', type=str, required=False, help='Select iteration number')
    args = parser.parse_args()

    return args

def plot_thermodynamic_mechanism(current_dir, protein, iteration_number, bins=40):
    """ Plot ordering versus reaction


    """
    # Get the Q files from the folding temperature directory                                                          
    path = current_dir + '/' + protein + '/Mut_'+iteration_number+'/'
    temp = open(path+'T_array_last.txt', 'r').readline().split('_')[0]
    # Tf_choice.txt should be manually edited to contain the temperature directory closest to the actual Tf 
    path += temp+'_1/'
    # We assume that each of the 3 directories for a given temperature 

    print "Loading Q.dat, Qres.dat"
    Q = np.loadtxt(path+"Q.dat")
    Qres = np.loadtxt(path+"qimap.dat")

    calculate_locals = True
    #This should soon be deprecated, since Qlocalres, etc. are no longer calculated
    try:
        Qres_local = np.loadtxt(path+"Qlocalres.dat")
    except: 
        calculate_locals = False

    if calculate_locals==True:
        Qres_nonlocal = np.loadtxt(path+"Qnonlocalres.dat")

    minQ = min(Q)
    maxQ = max(Q)
    incQ = (float(maxQ) - float(minQ))/bins

    Qprogress = np.zeros((bins,len(Qres[0])),float)
    Qprogress_local = np.zeros((bins,len(Qres[0])),float)
    Qprogress_nonlocal = np.zeros((bins,len(Qres[0])),float)
    counts = np.zeros(bins,float)
    counts_local = np.zeros(bins,float)
    counts_nonlocal = np.zeros(bins,float)
   
    print "Histogram Qres depending on Q"
    for i in range(len(Q)):
        for n in range(bins):
            if ((minQ + n*incQ)  < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                Qprogress[n,:] += Qres[i,:]
                counts[n] += 1
                if calculate_locals==True:
                    Qprogress_local[n,:] += Qres_local[i,:]
                    counts_local[n] += 1
                    Qprogress_nonlocal[n,:] += Qres_nonlocal[i,:]
                    counts_nonlocal[n] += 1

    native = Qres[0,:]
    native[ native == 0 ] = 1
    if calculate_locals ==True:
        native_local = Qres_local[0,:]
        native_local[ native_local == 0 ] = 1
        native_nonlocal = Qres_nonlocal[0,:]
        native_nonlocal[ native_nonlocal == 0 ] = 1

    print "Calculating fraction of formed native contacts for each bin"
    Qprogress = ((Qprogress/native).T)/counts
    if calculate_locals==True:
        Qprogress_local = ((Qprogress_local/native_local).T)/counts_local
        Qprogress_nonlocal = ((Qprogress_nonlocal/native_nonlocal).T)/counts_nonlocal

    print "Plotting thermodynamic mechanism"
    plt.figure()
    #plt.subplot(1,1,1,aspect=1)
    plt.pcolor(Qprogress,edgecolors='k')
    cbar = plt.colorbar()
    cbar.set_label("Fraction individual contacts formed $Q_{individual}$")
    plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
    plt.ylabel("Sequence index")
    plt.title("Thermodynamic Folding Progress")
    plt.savefig(current_dir+"/metrics/thermodynamic/thermodynamic_mechanism_profile_"+protein+"_"+iteration_number+".pdf")
    plt.clf()

    if calculate_locals==True:
        plt.figure()
        plt.pcolor(Qprogress_local,edgecolors='k')
        cbar = plt.colorbar()
        cbar.set_label("Fraction individual contacts formed $Q_{local}$")
        plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
        plt.ylabel("Sequence index")
        plt.title("Thermodynamic Folding Progress for $Q_{local}-protein "+protein)
        plt.savefig(current_dir+"/metrics/thermodynamic/thermodynamic_mechanism_profile_Qlocal_"+protein+
                    "-"+iteration_number+".pdf")
        plt.clf()

        plt.figure()
        plt.pcolor(Qprogress_nonlocal,edgecolors='k')
        cbar = plt.colorbar()
        cbar.set_label("Fraction individual contacts formed $Q_{nonlocal}$")
        plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
        plt.ylabel("Sequence index")
        plt.title("Thermodynamic Folding Progress for $Q_{nonlocal}$-protein "+protein)
        plt.savefig(current_dir+"/metrics/thermodynamic/thermodynamic_mechanism_profile_Qnonlocal_"+protein+
                    "_"+iteration_number+".pdf")
        plt.clf()


def get_beadbead_info(path):

    print " Loading BeadBead.dat"
    beadbead = np.loadtxt(path+"BeadBead.dat",dtype=str)
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int)
    pairs -= np.ones(pairs.shape,int)

    keep_interactions = np.zeros(len(interaction_numbers),int)
    for i in range(len(interaction_numbers)):
        if interaction_numbers[i] in ["ds","ss"]:
            pass
        else:
            keep_interactions[i] = int(interaction_numbers[i])

    sigij = sigij[keep_interactions != 0]
    epsij = epsij[keep_interactions != 0]
    deltaij = deltaij[keep_interactions != 0]
    pairs = pairs[keep_interactions != 0]

    return pairs, epsij, sigij


def get_contact_pair_distribution_versus_Q(current_dir, protein, iteration_number):
    # Get the BeadBead.dat file from the folding temperature directory
    path = current_dir + '/' + protein + '/Mut_'+iteration_number+'/'
    temp = open(path+'T_array_last.txt', 'r').readline().split('_')[0]
    path += temp+'_1/'

    pairs, epsij, sigij = get_beadbead_info(path)
    print " Loading trajectory"
    traj = md.load(path+"traj.xtc",top=path+'Native.pdb')
    print " Computing distances with mdtraj..."
    distances = md.compute_contacts(traj,pairs)
    print " Computing contacts with mdtraj..."
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)

    #Divide contacts into local and nonlocal
    local = []
    local_n = 0
    nonlocal_n = 0
    for i in range(len(pairs)):
        if abs(pairs[i][0]-pairs[i][1])<7:
            local.append(True)
            local_n+=1
        else:
            local.append(False)
            nonlocal_n +=1
    #Argghh had to do it the artisanal way
    local_contacts = np.zeros((len(contacts), local_n))
    nonlocal_contacts = np.zeros((len(contacts),nonlocal_n))
    local_n = 0
    nonlocal_n = 0

    for i in range(len(local)):
        if local[i]==True:
            local_contacts[:,local_n] = contacts[:,i]
            local_n+=1
        else:
            nonlocal_contacts[:,nonlocal_n]=contacts[:,i]
            nonlocal_n+=1
    #local_contacts = contacts[:,local==True]
    #nonlocal_contacts = contacts[:, local==False]

    print " Loading Q.dat"
    bins = 40
    Q = np.loadtxt(path+"Q.dat")
    Q /= float(max(Q))
    minQ = min(Q)
    maxQ = max(Q)
    incQ = (float(maxQ) - float(minQ))/bins
    
    # Added the local, nonlocal counters. 
    Q_i_for_bin = np.zeros((bins,len(contacts[0])),float)
    Qlocal_i_for_bin = np.zeros((bins,len(local_contacts[0])),float)
    Qnonlocal_i_for_bin = np.zeros((bins,len(nonlocal_contacts[0])),float)

    counts = np.zeros(bins,float)
    local_counts = np.zeros(bins, float)
    nonlocal_counts = np.zeros(bins, float)

    print " Histogram contacts depending on Q"
    for i in range(traj.n_frames):
        for n in range(bins):
            if ((minQ + n*incQ) < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                Q_i_for_bin[n,:] += contacts[i,:]
                Qlocal_i_for_bin[n,:] += local_contacts[i,:]
                Qnonlocal_i_for_bin[n,:] += nonlocal_contacts[i,:]
                counts[n] += 1
                local_counts[n] += 1
                nonlocal_counts[n] += 1

    print " Averaging contacts for each bin"
    Q_i_avg = (Q_i_for_bin.T/counts).T
    Qlocal_i_avg = (Qlocal_i_for_bin.T/local_counts).T
    Qnonlocal_i_avg = (Qnonlocal_i_for_bin.T/nonlocal_counts).T
    return minQ, maxQ, bins, Q_i_avg, Qlocal_i_avg, Qnonlocal_i_avg

def plot_contact_pair_distribution_versus_Q(current_dir, protein, iteration_number):

    minQ, maxQ, bins, Q_i_avg, Qlocal_i_avg, Qnonlocal_i_avg = get_contact_pair_distribution_versus_Q(current_dir, protein,
                                                                                                      iteration_number)
    
    # Same for all types of contacts

    Nbins = 25 
    Qbins = np.linspace(minQ,maxQ,bins-1)
    Qi_bins = np.linspace(0.,1.0,Nbins)
    Q_i_distribution = np.zeros((bins,Nbins-1))
    Qlocal_i_distribution = np.zeros((bins,Nbins-1))
    Qnonlocal_i_distribution = np.zeros((bins,Nbins-1))
    for j in range(0,bins):
        #n,bins,patches = plt.hist(Q_i_avg[j],bins=bins2,alpha=0.3,label=str(j),histtype='stepfilled')
        #plt.clf()
        n,tempbins = np.histogram(Q_i_avg[j],bins=Qi_bins,density=True)
        Q_i_distribution[j,:] = n
        n_local, tempbins_local = np.histogram(Qlocal_i_avg[j],bins=Qi_bins,density=True)
        Qlocal_i_distribution[j,:] = n_local
        n_nonlocal, tempbins_nonlocal = np.histogram(Qnonlocal_i_avg[j],bins=Qi_bins,density=True)
        Qnonlocal_i_distribution[j,:] = n_nonlocal
  
    x = np.linspace(minQ,maxQ,bins+1)
    y = np.linspace(0.,1.0,Nbins)
    X,Y = np.meshgrid(x,y)

    ## Set 
    Qi_dist = Q_i_distribution.T
    maxQi_dist = max(Qi_dist[1:-1,:].ravel())
    minQi_dist = min(Qi_dist[1:-1,:].ravel())

    Qlocal_i_dist = Qlocal_i_distribution.T
    maxQlocal_i_dist = max(Qlocal_i_dist[1:-1,:].ravel())
    minQlocal_i_dist = min(Qlocal_i_dist[1:-1,:].ravel())

    Qnonlocal_i_dist = Qnonlocal_i_distribution.T
    maxQnonlocal_i_dist = max(Qnonlocal_i_dist[1:-1,:].ravel())
    minQnonlocal_i_dist = min(Qnonlocal_i_dist[1:-1,:].ravel())
  
    print "  Saving data..."
    np.savetxt(current_dir+'/metrics/mechanism/Qi_distribution_'+protein+'_'+iteration_number+'.dat',Qi_dist)
    np.savetxt(current_dir+'/metrics/mechanism/Q_bins_'+protein+'_'+iteration_number+'.dat',x)
    np.savetxt(current_dir+'/metrics/mechanism/Qi_bins_'+protein+'_'+iteration_number+'.dat',y)
    np.savetxt(current_dir+'/metrics/mechanism/Qlocal_i_distribution_'+protein+'_'+iteration_number+'.dat', Qlocal_i_dist)
    np.savetxt(current_dir+'/metrics/mechanism/Qnonlocal_i_distribution_'+protein+'_'+iteration_number+'.dat', Qnonlocal_i_dist)

    print "  Plotting.."
    plt.figure()
    #plt.pcolor(X,Y,Q_i_distribution.T)
    #plt.pcolor(X,Y,Q_i_distribution.T[1:-1,:])
    plt.pcolor(X,Y,Qi_dist)
    plt.clim(minQi_dist,maxQi_dist)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Q_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity - Q -"+protein)
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    print "  Saving plot..."
    plt.savefig(current_dir+'/metrics/mechanism/mechanism_heterogeneity_versus_Q_'+protein+'_'+iteration_number+'.pdf')
    plt.clf()

    #Plot for Qlocal
    plt.figure() 
    plt.pcolor(X,Y,Qlocal_i_dist)
    plt.clim(minQlocal_i_dist,maxQlocal_i_dist)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Qlocal_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity - Qlocal -"+protein)
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    print "  Saving plot..."
    plt.savefig(current_dir+'/metrics/mechanism/mechanism_heterogeneity_local_versus_Q_'+protein+'_'+iteration_number+'.pdf')
    plt.clf()

    #Plot for Qnonlocal                                                                                                             
    plt.figure()
    plt.pcolor(X,Y,Qnonlocal_i_dist)
    plt.clim(minQnonlocal_i_dist,maxQnonlocal_i_dist)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Qnonlocal_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity - Qnonlocal -"+protein)
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    print "  Saving plot..."
    plt.savefig(current_dir+'/metrics/mechanism/mechanism_heterogeneity_nonlocal_versus_Q_'+protein+'_'+iteration_number+'.pdf')
    plt.clf()
    

if __name__ == "__main__":
    
    #List of possible proteins to choose from                                                                                   
    proteins_list = ['r15', 'r16', 'r17','1SHG', '1RIS', '1TEN', '1K85','1E0G','1E41','sh3']
    proteins_choice = [x for x in get_args(proteins_list).prot]
    #Mutation number (defaults to 0)                                                                                           
    iteration_number = get_args(proteins_list).iter
    
    if iteration_number == None:
        iteration_number = '0'
    # Determine if the metrics directory exists and create it if not                                                           
    current_dir = os.getcwd()
    # We will create a separate "metrics/mechanism" directory to store all the data                                                     
    if os.path.isdir(current_dir+'/metrics')==False:
        os.mkdir('metrics')
    
    os.chdir(current_dir+'/metrics')
    if os.path.isdir(current_dir+'/metrics/mechanism')==False:
        os.mkdir('mechanism')
    if os.path.isdir(current_dir+'/metrics/thermodynamic')==False:
        os.mkdir('thermodynamic')
    os.chdir(current_dir)

    for protein in proteins_choice:
        plot_contact_pair_distribution_versus_Q(current_dir, protein, iteration_number)
        plot_thermodynamic_mechanism(current_dir, protein, iteration_number, bins=40)
