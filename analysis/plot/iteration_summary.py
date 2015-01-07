''' Plot a summary of MC2004 iteration


Description:


'''


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

from mpl_toolkits.axes_grid1 import make_axes_locatable


colors = [('white')] + [(cm.Blues(i)) for i in xrange(1,256)]
global new_map
new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

colors2 = [('gray')] + [(cm.Blues(i)) for i in xrange(1,256)]
global new_map2
new_map2 = matplotlib.colors.LinearSegmentedColormap.from_list('new_map2', colors2, N=256)

#colors3 = [('gray')] + [(cm.Jet(i)) for i in xrange(1,256)]
#global new_map3
#new_map3 = matplotlib.colors.LinearSegmentedColormap.from_list('new_map3', colors3, N=256)

def get_contact_probability(name,iteration,n_residues,contacts,state_label,state_bound):

    cwd = os.getcwd()
    sub =  "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)
    
    if not os.path.exists("cont_prob_%s.dat" % state_label):
        contact_probability = np.zeros(len(contacts),float) 
        n_frames = 0.
        temps = [ x.rstrip("\n") for x in open("long_temps_last", "r").readlines() ]
        for i in range(len(temps)):
            T = temps[i]
            Q = np.loadtxt("%s/Q.dat" % T)
            qimap = np.loadtxt("%s/qimap.dat" % T)

            state_indicator = ((Q > state_bound[0]).astype(int)*(Q < state_bound[1]).astype(int)).astype(bool)
            n_frames += float(sum(state_indicator.astype(int)))
            contact_probability += sum(qimap[(state_indicator == True),:])
        contact_probability /= n_frames
        np.savetxt("cont_prob_%s.dat" % state_label,contact_probability)
    else:
        contact_probability = np.loadtxt("cont_prob_%s.dat" % state_label)

    os.chdir(cwd)
    return contact_probability

def get_iteration_data(name,iteration):
    ''' Get summary data for iteration '''
    n_residues = len(open("%s/Native.pdb" % name,"r").readlines()) - 1

    cwd = os.getcwd()
    sub =  "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)

    Tuse = open("long_temps_last","r").readlines()[0].rstrip("\n")
    Tf = float(open("long_Tf","r").read().rstrip("\n"))

    params = np.loadtxt("%s/pairwise_params" % Tuse,dtype=float)
    contacts = params[:,:2].astype(int)
    sign = params[:,3].astype(int)
    epsilons = np.loadtxt("%s/model_params" % Tuse,dtype=float)
    epsilons[sign == 3] = -1.*epsilons[sign == 3]

    epsilon_map = np.zeros((n_residues,n_residues))
    for i in range(len(contacts)):
        epsilon_map[contacts[i,1]-1,contacts[i,0]-1] = epsilons[i] 

    ## TODO:
    ##  - Should use pmf from data actually used in calculation not wham pmf.


    for line in open("short_wham/free.config","r"):
        if line.startswith("startTF"):
            startTF = float(line.split()[1])
        if line.startswith("deltaTF"):
            deltaTF = float(line.split()[1])
        if line.startswith("ntempsF"):
            ntempsF = int(line.split()[1])

    tempsF = [ float("%.1f" % (startTF + i*deltaTF)) for i in range(ntempsF) ]
    for j in range(ntempsF):
        if (Tf - tempsF[j]) <= 0:
            temp = "%.1f" % tempsF[j]
            whamT = "free" + temp.split(".")[0] + temp.split(".")[1]
            if os.path.exists("short_wham/%s" % whamT):
                whamFree = np.loadtxt("short_wham/%s" % whamT)
                break
            temp = "%.1f" % (tempsF[j] - 0.1)
            whamT = "free" + temp.split(".")[0] + temp.split(".")[1]
            if os.path.exists("short_wham/%s" % whamT):
                whamFree = np.loadtxt("short_wham/%s" % whamT)
                break
            temp = "%.1f" % (tempsF[j] + 0.1)
            whamT = "free" + temp.split(".")[0] + temp.split(".")[1]
            if os.path.exists("short_wham/%s" % whamT):
                whamFree = np.loadtxt("short_wham/%s" % whamT)
                break

    state_labels = []
    state_bounds = []
    for line in open("state_bounds.txt","r"):
        state_labels.append(line.split()[0])
        state_bounds.append([int(line.split()[1]),int(line.split()[2])])

    ddgsim = np.loadtxt("newton/sim_feature.dat")
    ddgsim_err = np.loadtxt("newton/sim_feature_err.dat")
    ddgexp = np.loadtxt("newton/target_feature.dat")
    ddgexp_err = np.loadtxt("newton/target_feature_err.dat")

    ddGsim = np.zeros((len(ddgsim),2),float)
    ddGsim[:,0] = ddgsim
    ddGsim[:,1] = ddgsim_err
    ddGexp = np.zeros((len(ddgexp),2),float)
    ddGexp[:,0] = ddgexp
    ddGexp[:,1] = ddgexp_err

    n_contacts = len(contacts)

    os.chdir(cwd)

    return epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, whamFree, state_labels, state_bounds, ddGsim, ddGexp

def plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=False):
    ''' Get  '''

    if iteration == 0:
        mineps = 0.
        maxeps = 2.

    if os.path.exists("%s/epsilon_range" % name):
        temp = np.loadtxt("%s/epsilon_range" % name) 
        mineps = temp[0]
        maxeps = temp[1]
    else:
        mineps = min(epsilons)
        maxeps = max(epsilons)

    #plt.pcolor(epsilon_map,cmap=new_map2,vmin=mineps,vmax=maxeps)      ## Blues cmap
    plt.pcolor(epsilon_map,cmap="jet",vmin=mineps,vmax=maxeps)
    cbar = plt.colorbar()
    cbar.set_clim(mineps,maxeps)
    plt.xticks(range(0,n_residues+1,10))
    plt.yticks(range(0,n_residues+1,10))
    if individual:
        plt.title("%s params after iteration %d" % (name,iteration))
    plt.grid(True)
    plt.xlim(0,n_residues)
    plt.ylim(0,n_residues)

def plot_epsilon_histogram(name,iteration,epsilons,individual=False):
    plt.hist(epsilons,bins=20,histtype="stepfilled",facecolor=new_map(0.75),alpha=0.4)
    if individual:
        plt.title("%s iteration %d $\\overline{\\epsilon} = %.3f$    $\\sigma^2 = %.3f$ " % \
                (name,iteration,np.mean(epsilons),np.std(epsilons)))
    else:
        plt.title("$\\overline{\\epsilon} = %.3f$    $\\sigma^2 = %.3f$ " % \
                (np.mean(epsilons),np.std(epsilons)))
    plt.xlabel("$\\epsilon$")
    plt.ylabel("Frequency")
    plt.grid(True)

def plot_free_energy(name,iteration,n_contacts,Tf,whamFree,state_labels,state_bounds,individual=False):
    TS = ((whamFree[:,0] > state_bounds[1][0]).astype(int)*(whamFree[:,0] < state_bounds[1][1]).astype(int)).astype(bool)
    barrier_height = max(whamFree[TS,1])

    plt.plot(whamFree[:,0],whamFree[:,1],'b',lw=2)
    plt.xlabel("Q")
    plt.ylabel("F(Q)")
    plt.xlim(0,n_contacts)
    plt.ylim(0,barrier_height+1)
    plt.vlines(state_bounds[0],0,1.5,lw=2)
    plt.vlines(state_bounds[1],barrier_height-1,barrier_height+0.5,lw=2)
    plt.vlines(state_bounds[2],0,1.5,lw=2)
    if individual:
        plt.title("%s iteration %d %.1f" % (name,iteration,Tf))
    else:
        plt.title("%.1f" % Tf)

def plot_ddG_comparison(name,iteration,ddGsim,ddGexp,individual=False):
    N = len(ddGsim)
    plt.errorbar(ddGsim[N/2:,0],ddGexp[N/2:,0],xerr=ddGsim[N/2:,1],yerr=ddGexp[N/2:,1],marker='o',linestyle="none",color='r',label="$\\Delta\\Delta G^{\\circ}$")
    plt.errorbar(ddGsim[:N/2,0],ddGexp[:N/2,0],xerr=ddGsim[:N/2,1],yerr=ddGexp[:N/2,1],marker='o',linestyle="none",color='b',label="$\\Delta\\Delta G^{\\dagger}$")
    plt.plot(range(-1,8),range(-1,8),color='k',lw=1)
    plt.xlabel("$\\Delta\\Delta G_{sim}$ (kT)")
    plt.ylabel("$\\Delta\\Delta G_{exp}$ (kT)")
    if individual:
        plt.title("%s iteration %d" % (name,iteration))
    maxddg = max([max(ddGsim[:,0]),max(ddGexp[:,0])])
    minddg = min([min(ddGsim[:,0]),min(ddGexp[:,0])])
    plt.xlim(minddg-0.5,maxddg+1)
    plt.ylim(minddg-0.5,maxddg+1)
    plt.grid(True)
    lg = plt.legend(loc=4)
    lg.draw_frame(False)

def plot_contact_probability(name,iteration,n_residues,contacts,state_label,state_bound,contact_probability,individual=False):
    ''' '''

    C = np.zeros((n_residues,n_residues),float)
    for j in range(len(contacts)):
        C[contacts[j,1]-1,contacts[j,0]-1] = contact_probability[j]

    #plt.pcolor(C,cmap=plt.get_cmap("Blues"),vmin=0.,vmax=1.0)
    plt.pcolor(C,cmap=new_map2,vmin=0.,vmax=1.0)
    if individual:
        cbar = plt.colorbar()
        cbar.set_clim(0,1)
    #for i in range(len(contacts)):
    #    plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map(contact_probability[i]))
    plt.xlim(0,n_residues)
    plt.ylim(0,n_residues)
    plt.xticks(range(0,n_residues,10))
    plt.yticks(range(0,n_residues,10))
    plt.grid(True)
    if individual:
        plt.title("%s iteration %d. %s contact probability" % (name,iteration,state_label))
    else:
        plt.title("%s contact probability" % (state_label))

def plot_contact_probability_subplot(name,iteration,n_residues,contacts,state_labels,Contact_maps):


    plt.figure(figsize=(9,4.2))
    for X in range(len(state_labels)):
        ax = plt.subplot(1,len(state_labels),X+1,aspect=1.0)
    
        C = np.zeros((n_residues,n_residues),float)
        for j in range(len(contacts)):
            C[contacts[j,1]-1,contacts[j,0]-1] = Contact_maps[X][j]

        #im = ax.pcolor(C,cmap=plt.get_cmap("Blues"),vmin=0.,vmax=1.0)
        im = ax.pcolor(C,cmap=new_map2,vmin=0.,vmax=1.0)
        plt.title("%s" % state_labels[X])
        #for i in range(len(contacts)):
        #    plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map(Contact_maps[X][i]))
        plt.xlim(0,n_residues)
        plt.ylim(0,n_residues)
        plt.xticks(range(0,n_residues,10))
        plt.yticks(range(0,n_residues,10))
        plt.grid(True)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_clim(0,1)

    plt.suptitle("%s iteration %d" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/contact_prob_all.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/contact_prob_all.png" % (name,iteration))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, whamFree, state_labels, state_bounds, ddGsim, ddGexp = get_iteration_data(name,iteration)
    
    print "Plotting summary for %s iteration %d..." % (name,iteration)

    if not os.path.exists("%s/iteration_%d/summary" % (name,iteration)):
        os.mkdir("%s/iteration_%d/summary" % (name,iteration))

    Contact_maps = []
    for X in range(len(state_labels)):
        plt.figure()
        contact_probability = get_contact_probability(name,iteration,n_residues,contacts,state_labels[X],state_bounds[X])
        plot_contact_probability(name,iteration,n_residues,contacts,state_labels[X],state_bounds[X],contact_probability,individual=True)
        Contact_maps.append(contact_probability)
        print " Saving: %s/iteration_%d/contact_prob_%s.pdf          - %s contact probabilities" % (name,iteration,state_labels[X],state_labels[X])
        plt.savefig("%s/iteration_%d/summary/contact_prob_%s.pdf" % (name,iteration,state_labels[X]))
        plt.savefig("%s/iteration_%d/summary/contact_prob_%s.png" % (name,iteration,state_labels[X]))
        plt.close()

    print " Saving subplot: %s/iteration_%d/summary/contact_prob_all.pdf       - contact probabilities" % (name,iteration)
    plot_contact_probability_subplot(name,iteration,n_residues,contacts,state_labels,Contact_maps)
    
    print "  Saving: %s/iteration_%d/summary/current_epsilon_map.pdf    - epsilon map" % (name,iteration)
    plt.figure()
    plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=True)
    plt.savefig("%s/iteration_%d/summary/current_epsilon_map.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/current_epsilon_map.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/summary/current_epsilon_hist.pdf   - epsilon histogram" % (name,iteration)
    plt.figure()
    plot_epsilon_histogram(name,iteration,epsilons,individual=True)
    plt.savefig("%s/iteration_%d/summary/current_epsilon_hist.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/current_epsilon_hist.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/summary/FreeEnergy_Q.pdf            - free energy" % (name,iteration)
    plt.figure()
    plot_free_energy(name,iteration,n_contacts,Tf,whamFree,state_labels,state_bounds,individual=True)
    plt.savefig("%s/iteration_%d/summary/FreeEnergy_Q.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/FreeEnergy_Q.png" % (name,iteration))
    plt.close()

    print "  Saving: %s/iteration_%d/compareddG.pdf              - ddG comparison" % (name,iteration)
    plt.figure()
    plot_ddG_comparison(name,iteration,ddGsim,ddGexp,individual=True)
    plt.savefig("%s/iteration_%d/summary/compareddG.pdf" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.png" % (name,iteration,name,iteration))
    plt.close()

    print " Summary figure: %s/iteration_%d/summary_%s_%d.pdf" % (name,iteration,name,iteration)
    plt.figure(figsize=(12,10))
    plt.subplot(2,2,1)
    plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues)
    plt.subplot(2,2,2)
    plot_epsilon_histogram(name,iteration,epsilons)
    plt.subplot(2,2,3)
    plot_free_energy(name,iteration,n_contacts,Tf,whamFree,state_labels,state_bounds)
    plt.subplot(2,2,4)
    plot_ddG_comparison(name,iteration,ddGsim,ddGexp)
    plt.suptitle("%s iteration %d" % (name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/summary/summary_%s_%d.png" % (name,iteration,name,iteration))
    #plt.show()
