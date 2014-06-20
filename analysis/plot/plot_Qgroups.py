""" Different ways to visualize the folding mechanism

Description:

    There are many ways to visualize the folding mechanism. For a description
of several topological measures see reference (1).


To Do: 
    
- Have option to select particular partitions of contacts



References:

(1) Chavez, L.; Onuchic, J.; Clementi, C. Quantifying the roughness on the
free energy landscape: Entropic bottlenecks and protein folding rates. J.
Am. Chem. Soc. 2004, 126, 8426-8432.
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
    """ Plot residue ordering versus folding progress.

    Description:

        This plots the average native structure for each residue as a function
    of folding progress (overall Q). This visualizes how ordering proceeds at
    different parts of the chain.

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

def get_contact_groups_distribution_versus_Q():
    """ Calculate contact probabilities and histogram by foldedness (Q)

    Description:
        
        Return a contact probabilities


    """
    #groups = []
    #for line in open("","r"):
    #    pass    

    pairs, epsij, sigij = get_beadbead_info()
    print " Loading trajectory"
    traj = md.load("traj.xtc",top="Native.pdb")
    print " Computing distances with mdtraj..."
    distances = md.compute_contacts(traj,pairs)
    print " Computing contacts with mdtraj..."
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)


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

def get_contact_pair_distribution_versus_Q():
    """ Calculate contact probabilities and histogram by foldedness (Q)

    Description:
        
        Return a contact probabilities


    """
    pairs, epsij, sigij = get_beadbead_info()
    print " Loading trajectory"
    traj = md.load("traj.xtc",top="Native.pdb")
    print " Computing distances with mdtraj..."
    distances = md.compute_contacts(traj,pairs)
    print " Computing contacts with mdtraj..."
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)


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

def plot_contact_groups_versus_Q():
    """ Plot contact group ordering versus folding progress.

    Description:

    """

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

    ## Set 
    Qi_dist = Q_i_distribution.T
    maxQi_dist = max(Qi_dist[1:-1,:].ravel())
    minQi_dist = min(Qi_dist[1:-1,:].ravel())

    print "  Saving data..."
    np.savetxt("plots/Qi_distribution.dat",Qi_dist)
    np.savetxt("plots/Q_bins.dat",x)
    np.savetxt("plots/Qi_bins.dat",y)


    print "  Plotting.."
    plt.figure()
    #plt.pcolor(X,Y,Q_i_distribution.T)
    #plt.pcolor(X,Y,Q_i_distribution.T[1:-1,:])
    plt.pcolor(X,Y,Qi_dist)
    plt.clim(minQi_dist,maxQi_dist)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Q_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity")
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    print "  Saving plot..."
    plt.savefig("plots/mechanism_heterogeneity_versus_Q.pdf")
    plt.show()

def plot_contact_pair_distribution_versus_Q():
    """ Plot contact ordering versus folding progress.

    Description:

        This plots the distribution of contact probabilities as a function
    of folding progress (overall Q). This is useful for discerning if ordering
    is uniform, i.e. all residues order gradually and together, or 
    heterogeneous, i.e. a set of contacts form a higher probability. 
        The distribution of contact probabilities for a given overall Q has a 
    mean of the overall Q. Therefore a uniform mechanism has a unimodal 
    distribution whose mean follows the y=x line. Whereas a hetergeneous 
    mechanism may have a bimodal (or very spread out) distribution where one
    band of contacts form early before/during transition state and one band of
    contacts form after.
        The variance of a slice at a given Q divided by the maximum possible 
    variance Q*(1-Q) is the route measure at that Q, R(Q). See reference (1)
    on route measure.
        This is related to the route entropy discussed by Plotkin et.al. in 
    reference (2). This plot is like ref (2) fig 2, but with all possible
    contacts plotted and taking as a density.


    References:
    (1) Chavez, L.; Onuchic, J.; Clementi, C. Quantifying the roughness on the
    free energy landscape: Entropic bottlenecks and protein folding rates. J.
    Am. Chem. Soc. 2004, 126, 8426-8432.

    (2) Plotkin, S.; Onuchic, J. Structural and energetic heterogeneity 
    in protein folding. I. Theory. J. Chem. Phys. 2002, 116, 5263.
    """

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

    ## Set 
    Qi_dist = Q_i_distribution.T
    maxQi_dist = max(Qi_dist[1:-1,:].ravel())
    minQi_dist = min(Qi_dist[1:-1,:].ravel())

    print "  Saving data..."
    np.savetxt("plots/Qi_distribution.dat",Qi_dist)
    np.savetxt("plots/Q_bins.dat",x)
    np.savetxt("plots/Qi_bins.dat",y)


    print "  Plotting.."
    plt.figure()
    #plt.pcolor(X,Y,Q_i_distribution.T)
    #plt.pcolor(X,Y,Q_i_distribution.T[1:-1,:])
    plt.pcolor(X,Y,Qi_dist)
    plt.clim(minQi_dist,maxQi_dist)
    plt.xlabel("Folding Reaction ($Q$)")
    plt.ylabel("$\\left< Q_i \\right>$ Distribution")
    plt.title("Mechanism Heterogeneity")
    plt.xlim(minQ,maxQ)
    plt.ylim(0.,1)
    print "  Saving plot..."
    plt.savefig("plots/mechanism_heterogeneity_versus_Q.pdf")
    plt.show()


def plot_route_measure_versus_Q():
    """ Plot route measure as defined in ref (1)

    Description:

        Given the distribution of contact probabilities at a given overall
    degree of foldedness (Q), the route measure is the variance of this 
    distribution divided by the maximum possible variance at the given Q, which
    is Q*(1-Q).
        A route measure R(Q)=1 means that contacts are forming in a binary 
    manner, i.e. contact probabilities take only values 0 or 1. On the other
    hand R(Q)=0 means every contact probabilities is equal to Q, i.e. the 
    mechanism is uniform.
        See reference (1) for more information.


    References:
    (1) Chavez, L.; Onuchic, J.; Clementi, C. Quantifying the roughness on the
    free energy landscape: Entropic bottlenecks and protein folding rates. J.
    Am. Chem. Soc. 2004, 126, 8426-8432.

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

if __name__ == "__main__":
    if not os.path.exists("plots"):
        os.mkdir("plots")
    #plot_thermodynamic_mechanism(5)
    #plot_contact_pair_distribution_versus_Q()
    #plot_route_measure_versus_Q()
    #plot_contact_pair_distribution_and_route_measure_versus_Q()

    element_bounds = [0,5,13,21,32,41,49,56]
    elements = [ "beta_"+str(i+1) for i in range(len(element_bounds)-1) ]
    Qgroups = {}
    Qref = np.loadtxt("../../Qref_cryst.dat")

    for i in range(len(Qref)):
        for j in range(len(Qref)):
            if Qref[i,j] == 1:
                elm1 = -1
                elm2 = -1
                for k in range(len(element_bounds)-1):
                    if (i >= element_bounds[k]) and (i < element_bounds[k+1]):
                        elm1 = k
                    if (j >= element_bounds[k]) and (j < element_bounds[k+1]):
                        elm2 = k
                if (elm1 == -1) or (elm2 == -1):
                    pass
                else:
                    group = elements[elm1]+"-"+elements[elm2]
                    if Qgroups.has_key(group):
                        Qgroups[group].append([i,j])
                    else:
                        Qgroups[group] = [[i,j]]

    keys = Qgroups.keys()
    for key in keys:
        if len(Qgroups[key]) <= 5:
            small_group = Qgroups.pop(key,None)
            continue

            #for pair in small_group:
            #    min_dist = 100000
            #    new_keys = Qgroups.keys()
            #    for new_key in new_keys: 
            #        for pair2 in Qgroups[new_key]:
            #            dist = np.sqrt((pair[0]-pair2[0])**2 + (pair[1]-pair2[1])**2)
            #            if dist < min_dist:
            #                future_group = new_key
            #    Qgroups[future_group].append(pair)

                

    Qgroups_string = ""
    for key in Qgroups.keys():
        Qgroups_string += "[ "+key+"]\n"
        for contact in Qgroups[key]:
            Qgroups_string += "%d %d\n" % (contact[0],contact[1])

    open("Qgroups.dat","w").write(Qgroups_string)
        

    colors = ["blue","green","red","indigo","maroon","chocolate","salmon","darkviolet","lightgreen","sienna","violet","orchid"]

    for i in range(len(Qgroups.keys())):
        key = Qgroups.keys()[i]
        j = 0
        label = "$\\"+key.split("-")[0]+"\\"+key.split("-")[1]+"$"
        for pair in Qgroups[key]:
            if j == 0:
                j += 1
                plt.plot(pair[1]+1,pair[0]+1,marker='s',markersize=6,markeredgecolor=colors[i],color=colors[i],label=label)
            else:
                plt.plot(pair[1]+1,pair[0]+1,marker='s',markersize=6,markeredgecolor=colors[i],color=colors[i])

    plt.legend(loc=2)
    plt.xlim(0,len(Qref))
    plt.ylim(0,len(Qref))
    plt.show()


