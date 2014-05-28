""" Different ways to visualize the folding mechanism

Description:

    There are many ways to visualize the folding mechanism.

"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

import mdtraj as md

def plot_kinetic_mechanism():
    ''' The kinetic mechanism is defined as the mechanism of only successful or
        "reactive" trajectories (i.e. those that traverse from unfolded to
        folded before returning to unfolded). This trims away unreactive 
        fluctuactions from each of the states.

        NOT DONE.'''
    coord = "Q"
    statefile = open(coord+"_states.txt","r").readlines()[1:]
    bounds = []
    for line in statefile:
        bounds.append([line.split()[0],float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    print "Loading Q.dat, Qres.dat"
    Q = np.loadtxt("Q.dat")
    Qres = np.loadtxt("Qres.dat")


    minQ = min(Q)
    maxQ = max(Q)
    bins = 50
    incQ = (float(maxQ) - float(minQ))/bins

    Qprogress = np.zeros((bins,len(Qres[0])),float)
    counts = np.zeros(bins,float)

    left_bound = bounds[0][3]
    right_bound = bounds[2][2]
    folding = 0
    for i in range(1,len(Q)):

        if folding == 0:
            if (Q[i] > left_bound) and (Q[i-1] < left_bound):
                folding = 1

def plot_thermodynamic_mechanism(bins=10):
    """ Plot ordering versus reaction


    """
    bins = 10

    print "Loading Q.dat, Qres.dat"
    Q = np.loadtxt("Q.dat")
    Qres = np.loadtxt("Qres.dat")

    minQ = min(Q)
    maxQ = max(Q)
    incQ = (float(maxQ) - float(minQ))/bins

    Qprogress = np.zeros((bins,len(Qres[0])),float)
    counts = np.zeros(bins,float)


    print "Histogram Qres depending on Q"
    for i in range(len(Q)):
        for n in range(bins):
            if ((minQ + n*incQ)  < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                Qprogress[n,:] += Qres[i,:]
                counts[n] += 1

    native = Qres[0,:]
    native[ native == 0 ] = 1

    print "Calculating fraction of formed native contacts for each bin"
    Qprogress = ((Qprogress/native).T)/counts

    print "Plotting thermodynamic mechanism"
    plt.figure()
    #plt.subplot(1,1,1,aspect=1)
    plt.pcolor(Qprogress,edgecolors='k')
    cbar = plt.colorbar()
    cbar.set_label("Fraction local contacts formed $Q_{local}$")
    plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
    plt.ylabel("Sequence index")
    plt.title("Thermodynamic Folding Progress")
    plt.savefig("thermodynamic_mechanism_profile.pdf")
    #plt.savefig("mechanism_profile_square.pdf")
    plt.show()


def get_beadbead_info():

    print " Loading BeadBead.dat"
    beadbead = np.loadtxt("BeadBead.dat",dtype=str)
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


def get_contact_pair_distribution_versus_Q():
    pairs, epsij, sigij = get_beadbead_info()
    print " Loading trajectory"
    traj = md.load("traj.xtc",top="Native.pdb")
    print " Computing distances with mdtraj..."
    distances = md.compute_contacts(traj,pairs)
    print " Computing contacts with mdtraj..."
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)

    if not os.path.exists("plots"):
        os.mkdir("plots")

    print " Loading Q.dat"
    bins = 40
    Q = np.loadtxt("Q.dat")
    Q /= float(max(Q))
    minQ = min(Q)
    maxQ = max(Q)
    incQ = (float(maxQ) - float(minQ))/bins

    Q_i_for_bin = np.zeros((bins,len(contacts[0])),float)
    counts = np.zeros(bins,float)

    print " Histogram contacts depending on Q"
    for i in range(traj.n_frames):
        for n in range(bins):
            if ((minQ + n*incQ) < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                Q_i_for_bin[n,:] += contacts[i,:]
                counts[n] += 1

    print " Averaging contacts for each bin"
    Q_i_avg = (Q_i_for_bin.T/counts).T
    return minQ, maxQ, bins, Q_i_avg

def plot_contact_pair_distribution_versus_Q():

    minQ, maxQ, bins, Q_i_avg = get_contact_pair_distribution_versus_Q()

    Nbins = 25 
    Qbins = np.linspace(minQ,maxQ,bins-1)
    Qi_bins = np.linspace(0.,1.0,Nbins)
    Q_i_distribution = np.zeros((bins,Nbins-1))
    for j in range(0,bins):
        #n,bins,patches = plt.hist(Q_i_avg[j],bins=bins2,alpha=0.3,label=str(j),histtype='stepfilled')
        #plt.clf()
        n,tempbins = np.histogram(Q_i_avg[j],bins=Qi_bins,density=True)
        Q_i_distribution[j,:] = n

    x = np.linspace(minQ,maxQ,bins+1)
    y = np.linspace(0.,1.0,Nbins)
    X,Y = np.meshgrid(x,y)

    print " Plotting.."
    plt.figure()
    plt.pcolor(X,Y,Q_i_distribution.T)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Q_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity")
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    plt.savefig("plots/mechanism_heterogeneity_versus_Q.pdf")
    plt.show()


def plot_route_measure_versus_Q():
    """ Plot 

    """

    minQ, maxQ, N_Qbins, Q_i_avg = get_contact_pair_distribution_versus_Q()

    print " Computing route measure"
    incQ = (float(maxQ) - float(minQ))/N_Qbins
    route_measure = np.zeros(N_Qbins,float)
    for k in range(N_Qbins):
        qtemp = minQ + k*incQ
        R_Q = (1./(qtemp*(1.-qtemp)))*np.mean((Q_i_avg[k] - qtemp)**2)
        route_measure[k] = R_Q

    print " Plotting.."
    Qbins = np.arange(N_Qbins)*incQ + minQ
    plt.plot(Qbins,route_measure,'ro')
    plt.title("Route Measure")
    plt.xlabel("Q")
    plt.ylabel("R(Q)")
    plt.xlim(0,1)
    plt.ylim(0,0.5)
    #plt.ylim(0,max(route_measure)+0.05)
    print " Saving.."
    plt.savefig("plots/route_measure.pdf")
    plt.show()

def plot_contact_pair_distribution_and_route_measure_versus_Q():

    minQ, maxQ, N_Qbins, Q_i_avg = get_contact_pair_distribution_versus_Q()

    N_Qi_bins = 25 
    Qbins = np.linspace(minQ,maxQ,N_Qbins-1)
    Qi_bins = np.linspace(0.,1.0,N_Qi_bins)
    Q_i_distribution = np.zeros((N_Qbins,N_Qi_bins-1))
    for j in range(0,N_Qbins):
        #n,bins,patches = plt.hist(Q_i_avg[j],bins=bins2,alpha=0.3,label=str(j),histtype='stepfilled')
        #plt.clf()
        n,tempbins = np.histogram(Q_i_avg[j],bins=Qi_bins,density=True)
        Q_i_distribution[j,:] = n

    x = np.linspace(minQ,maxQ,N_Qbins+1)
    y = np.linspace(0.,1.0,N_Qi_bins)
    X,Y = np.meshgrid(x,y)

    print " Computing route measure"
    incQ = (float(maxQ) - float(minQ))/N_Qbins
    route_measure = np.zeros(bins,float)
    for k in range(bins):
        qtemp = minQ + k*incQ
        R_Q = (1./(qtemp*(1.-qtemp)))*np.mean((Q_i_avg[k] - qtemp)**2)
        route_measure[k] = R_Q

    print " Plotting.."
    #Qbins = np.arange(N_Qbins)*incQ + minQ

if __name__ == "__main__":
    pass
    #plot_thermodynamic_mechanism(5)
    #plot_contact_pair_distribution_versus_Q()
    plot_route_measure_versus_Q()
    #plot_contact_pair_distribution_and_route_measure_versus_Q()

