''' Calcuate and plot contact probability (Qi) vs foldedness (Q) 

Description:
    Contact probabilities are the elementary descriptors of the folding
mechanism. Since there are usually >100 contacts it is usually too hard to
visualize/understand whats going on just by looking at the contact
probabilities. 
    This module first calculates then plots several coarse ways of visualizing
contact probabilities:

1. Route measure (see reference (1))
    The normalized variance of the contact probilities as a function of Q (the
mean of the contact probailities at a given Q is equal to Q). The variance
indicates how polarized the mechanism is: a route measure of 1 is the maximum
possible polarization where contacts are formed in a binary way (1 or 0), a
lower route measure indicates that contact probabilities are fractional due to
averaging of diverse structures.

2. Qi hist vs Q
    More details can be seen by visualizing how the distribution of contact
probabilities change as a function of Q. For example, is it bimodial with
clusters of well-formed and un-formed contacts? or is it unimodal?
Hypothetically different distributions can result in the same route measure,
e.g. a vary wide unimodal distribution and a tightly clustered bimodal
distribution. 
    For this figure, contact probabilities are histogrammed at a given Q
then the values of the histogram are visualized in a color plot.

3. Qi vs Q colored by loop length
    The greatest level of detail is to look at all individual contact
probabilities as a function of Q. This is not easy to digest visually so this
figure separates loop lengths (less)greater than the average loop length into
two adjacent subplots (contact loop length is the sequence separation of
residues in contact). The contact probabilities are also colored by loop length
from shortest loops with shortest wavelength (black/blue) to longest loops
with the longest wavelength (red/white).


See also:
1. Qgroups_ss
2. Qgroups
3. mechanism

References:

(1) Chavez, L.; Onuchic, J.; Clementi, C. Quantifying the roughness on the
free energy landscape: Entropic bottlenecks and protein folding rates. J.
Am. Chem. Soc. 2004, 126, 8426-8432.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import argparse

import extract_params 

global SKIP_INTERACTIONS
SKIP_INTERACTIONS = [1,8,9]


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

def plot_QivsQ(name,iteration,Qbins,Qi_vs_Q,n_bins,epsilons,loops,state_bounds):
    ''' Plot contact ordering versus folding progress.

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
    '''

    ## Contact probability colored by loop length
    sortedindx = loops.argsort()
    fig, axes = plt.subplots(1,2,sharey=True,figsize=(11,5.5))
    for i in range(len(loops)):
        index = sortedindx[i]
        if (loops[index] > np.mean(loops)):
            ax = axes[0]
        else:
            ax = axes[1]
        ax.plot(Qbins,Qi_vs_Q[:,index],lw=2.5,alpha=0.5,color=cm.spectral((loops[index] - min(loops))/float(max(loops))))

    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("Contact Probability $\\langle Q_i \\rangle$")
    axes[0].set_title("$l_i > \\overline{l}$")
    axes[1].set_title("$l_i < \\overline{l}$")
    fig.text(0.5,0.04,"Foldedness $Q$",ha='center',va='center')
    fig.text(0.5,0.96,"Longer loops = longer wavelengths. %s iteration %d" % (name,iteration),ha='center',va='center')
    plt.savefig("%s/iteration_%d/plots/QivsQ_loops_%s_%d.png" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QivsQ_loops_%s_%d.pdf" % (name,iteration,name,iteration))

    ## Contact probability growth colored by loop length
    dQi_dQ = (np.gradient(Qi_vs_Q,Qbins[1]-Qbins[0],Qbins[1]-Qbins[0]))[0]
    fig, axes = plt.subplots(1,2,sharey=True,figsize=(11,5.5))
    for i in range(len(loops)):
        index = sortedindx[i]
        if (loops[index] > np.mean(loops)):
            ax = axes[0]
        else:
            ax = axes[1]
        ax.plot(Qbins,dQi_dQ[:,index],lw=2.5,alpha=0.5,color=cm.spectral((loops[index] - min(loops))/float(max(loops))))

    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("Contact Growth $\\frac{d\\langle Q_i \\rangle}{dQ}$")
    axes[0].set_title("$l_i > \\overline{l}$")
    axes[1].set_title("$l_i < \\overline{l}$")
    fig.text(0.5,0.04,"Foldedness $Q$",ha='center',va='center')
    fig.text(0.5,0.96,"Longer loops = longer wavelengths. %s iteration %d" % (name,iteration),ha='center',va='center')
    plt.savefig("%s/iteration_%d/plots/dQidQ_loops_%s_%d.png" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/dQidQ_loops_%s_%d.pdf" % (name,iteration,name,iteration))

def some_old_plot_Qgroups():
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

def get_some_iteration_data(name,iteration,n_bins):
    ''' Get summary data for iteration '''
    Tuse = open("%s/iteration_%d/long_temps_last" % (name,iteration),"r").readlines()[0].rstrip("\n")


    pair_file = "%s/iteration_%d/%s/pairwise_params" % (name,iteration,Tuse)
    param_file = "%s/iteration_%d/%s/model_params" % (name,iteration,Tuse)

    contacts, epsilons = extract_params.get_contacts_epsilons(pair_file,param_file)

    loops = contacts[:,1] - contacts[:,0]
    n_residues = len(open("%s/Native.pdb" % name,"r").readlines()) - 1

    state_labels = []
    state_bounds = []
    for line in open("%s/iteration_%d/state_bounds.txt" % (name,iteration),"r"):
        state_labels.append(line.split()[0])
        state_bounds.append([int(line.split()[1]),int(line.split()[2])])

    n_contacts = len(contacts)

    ## Calculate or load cont. prob. Qi versus reaction coordinate, Q.
    Qbins, Qi_vs_Q = get_contact_probability_versus_Q(name,iteration,n_bins)

    return epsilons, loops, n_residues, contacts, n_contacts, state_labels, state_bounds, Qbins, Qi_vs_Q

def get_contact_probability_versus_Q(name,iteration,n_bins):
    ''' Get the contact probabilities versus Q

    Description:
        
        If the contact probabilities have not been calculated for this
    name and iteration, QivsQ.dat, then calculate them and save them.
    
    '''
    print "Plotting mechanism info for %s iteration %d..." % (name,iteration)
    if not os.path.exists("%s/iteration_%d/QivsQ.dat" % (name,iteration)):
        print "  calculating Qi vs Q"
        n_frames = 0.
        temps = [ x.rstrip("\n") for x in open("%s/iteration_%d/long_temps_last" % (name,iteration), "r").readlines() ]
        for i in range(len(temps)):
            T = temps[i]
            Q_temp = np.loadtxt("%s/iteration_%d/%s/Q.dat" % (name,iteration,T))
            Qi_temp = np.loadtxt("%s/iteration_%d/%s/qimap.dat" % (name,iteration,T))
            if i == 0:
                Qi = Qi_temp
                Q = Q_temp
            else:
                Qi = np.concatenate((Qi,Qi_temp),axis=0)
                Q = np.concatenate((Q,Q_temp),axis=0)

        counts = np.zeros(n_bins)
        Qi_vs_bins = np.zeros((n_bins,len(Qi[0,:])),float)
        minQ = min(Q)
        maxQ = max(Q)
        incQ = (float(maxQ) - float(minQ))/float(n_bins)

        print "  sorting Qi"
        for i in range(len(Q)):
            for n in range(n_bins):
                if ((minQ + n*incQ) < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                    Qi_vs_bins[n,:] += Qi[i,:]
                    counts[n] += 1.

        Qi_vs_Q = (Qi_vs_bins.T/counts).T
        Qbins = np.linspace(minQ,maxQ,n_bins)
        np.savetxt("%s/iteration_%d/QivsQ.dat" % (name,iteration),Qi_vs_Q)
        np.savetxt("%s/iteration_%d/Qbins.dat" % (name,iteration),Qbins)

    else:
        print "  loading Qi vs Q"
        Qi_vs_Q = np.loadtxt("%s/iteration_%d/QivsQ.dat" % (name,iteration))
        Qbins = np.loadtxt("%s/iteration_%d/Qbins.dat" % (name,iteration))
    return Qbins, Qi_vs_Q

def plot_route_measure(name,iteration,Qbins,Qi_vs_Q,n_bins):
    ''' Plot route measure as defined in ref (1)

    Description:

        Given the distribution of contact probabilities at a given overall
    degree of foldedness (Q), the route measure is the variance of this 
    distribution divided by the maximum possible variance at the given Q, which
    is Q*(1-Q).
        A route measure R(Q)=1 means that contacts are forming in a binary 
    manner, i.e. contact probabilities take only values 0 or 1. On the other
    hand R(Q)=0 means every contact probabilities is equal to Q, i.e. the 
    mechanism is uniform.
        See reference (1) for more information. And reference (2) for original
    discussion.


    References:
    (1) Chavez, L.; Onuchic, J.; Clementi, C. Quantifying the roughness on the
    free energy landscape: Entropic bottlenecks and protein folding rates. J.
    Am. Chem. Soc. 2004, 126, 8426-8432.

    (2) Plotkin, S.; Onuchic, J.; Structural and energetic heterogeneity in 
    protein folding. I, Theory. Jour. Chem. Phys. 2002, 116, 12. 5263-5283

    '''

    Q = Qbins/float(max(Qbins))
    route = np.zeros(n_bins)
    for i in range(n_bins):
        if (Q[i] == 0) or (Q[i] == 1):
            pass
        else:
            route[i] = (1./(Q[i]*(1. - Q[i])))*(np.std(Qi_vs_Q[i,:])**2)

    np.savetxt("%s/iteration_%d/route.dat" % (name,iteration),route)
    plt.figure()
    plt.plot(Q,route,lw=2,color='r')
    plt.xlabel("Q")
    plt.ylabel("R(Q)")
    plt.title("Route measure %s iteration %d"  % (name,iteration))
    plt.savefig("%s/iteration_%d/plots/route_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/route_%s_%d.png" % (name,iteration,name,iteration))

def plot_Qi_histogram_vs_Q(name,iteration,Qbins,Qi_vs_Q,n_bins):

    n_probbins = 20
    probbins = np.linspace(0,1,n_probbins+1)
    deltaprob = probbins[1] - probbins[0]
    Qihist_vs_Q = np.zeros((n_probbins,n_bins),float)
    for i in range(n_bins):
        n, bins = np.histogram(Qi_vs_Q[i,:],density=True,bins=probbins)
        Qihist_vs_Q[:,i] = n*deltaprob

    Q = Qbins/float(max(Qbins))
    X,Y = np.meshgrid(Q,probbins)
    minQprob = min(Qihist_vs_Q[1:-1,:].ravel())
    maxQprob = max(Qihist_vs_Q[1:-1,:].ravel())

    plt.figure()
    #plt.pcolor(X[1:-1,:],Y[1:-1,:],Qihist_vs_Q[1:-1,:])
    plt.pcolor(X,Y,Qihist_vs_Q,vmin=minQprob,vmax=maxQprob)
    cbar = plt.colorbar()
    #cbar.set_clim(0.2*deltaprob,0.8*deltaprob)
    cbar.set_clim(minQprob,maxQprob)
    cbar.set_label("Probability")
    plt.ylim(0,max(Y.ravel()))
    plt.xlim(0,max(X.ravel()))
    plt.xlabel("Folding Reaction $Q$")
    plt.ylabel("$\\left< Q_i \\right>$ Distribution")
    plt.title("Contact formation %s iteration %s"  % (name,iteration))
    plt.savefig("%s/iteration_%d/plots/QihistvsQ_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QihistvsQ_%s_%d.png" % (name,iteration,name,iteration))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    parser.add_argument('--n_bins', type=int, default=35, help='Number of bins along Q.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    n_bins = args.n_bins

    if not os.path.exists("%s/iteration_%d/plots" % (name,iteration)):
        os.mkdir("%s/iteration_%d/plots" % (name,iteration))

    ## Get some iteration data
    epsilons, loops, n_residues, contacts, n_contacts, state_labels, state_bounds, Qbins, Qi_vs_Q = get_some_iteration_data(name,iteration,n_bins)

    print "  plotting route measure.     saving %s/iteration_%d/plots/route_%s_%d.png" % (name,iteration,name,iteration)
    plot_route_measure(name,iteration,Qbins,Qi_vs_Q,n_bins)

    print "  plotting Qi histogram vs Q. saving %s/iteration_%d/plots/QihistvsQ_%s_%d.png" % (name,iteration,name,iteration)
    plot_Qi_histogram_vs_Q(name,iteration,Qbins,Qi_vs_Q,n_bins)

    print "  plotting Qi vs Q            saving %s/iteration_%d/plots/QivsQ_loops_%s_%d.png" % (name,iteration,name,iteration)
    print "  plotting dQi/dQ vs Q        saving %s/iteration_%d/plots/dQidQ_loops_%s_%d.png" % (name,iteration,name,iteration)
    plot_QivsQ(name,iteration,Qbins,Qi_vs_Q,n_bins,epsilons,loops,state_bounds)

    #plt.show()


