import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib

import argparse

colors = [('white')] + [(cm.Blues(i)) for i in xrange(1,256)]
global new_map
new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

def get_iteration_data(name,iteration):
    """ Get summary data for iteration """
    Tuse = open("%s/Mut_%d/T_array_last.txt" % (name,iteration),"r").readlines()[0].rstrip("\n")
    Tf = float(open("%s/Mut_%d/Tf.txt" % (name,iteration),"r").read().rstrip("\n"))

    beadbead = np.loadtxt("%s/Mut_%d/%s/BeadBead.dat" % (name,iteration,Tuse),usecols=(0,1,6))
    contacts = beadbead[:,:2]
    epsilons = beadbead[:,2]

    n_residues = len(open("%s/Native.pdb" % name,"r").readlines()) - 1
    epsilon_map = np.zeros((n_residues,n_residues))
    for i in range(len(contacts)):
        epsilon_map[contacts[i,1]-1,contacts[i,0]-1] = epsilons[i] 

    for line in open("%s/Mut_%d/whamQ/free.config" % (name,iteration),"r"):
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
            if os.path.exists("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT)):
                whamFree = np.loadtxt("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT))
                break
            temp = "%.1f" % (tempsF[j] - 0.1)
            whamT = "free" + temp.split(".")[0] + temp.split(".")[1]
            if os.path.exists("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT)):
                whamFree = np.loadtxt("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT))
                break
            temp = "%.1f" % (tempsF[j] + 0.1)
            whamT = "free" + temp.split(".")[0] + temp.split(".")[1]
            if os.path.exists("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT)):
                whamFree = np.loadtxt("%s/Mut_%d/whamQ/%s" % (name,iteration,whamT))
                break

    state_labels = []
    state_bounds = []
    for line in open("%s/Mut_%d/state_bounds.txt" % (name,iteration),"r"):
        state_labels.append(line.split()[0])
        state_bounds.append([int(line.split()[1]),int(line.split()[2])])

    ddgsim = np.loadtxt("%s/Mut_%d/mut/ddGsim.dat" % (name,iteration))
    ddgsim_err = np.loadtxt("%s/Mut_%d/mut/ddGsim_err.dat" % (name,iteration))
    ddgexp = np.loadtxt("%s/Mut_%d/mut/ddGexp.dat" % (name,iteration))
    ddgexp_err = np.loadtxt("%s/Mut_%d/mut/ddGexp_err.dat" % (name,iteration))

    ddGsim = np.zeros((len(ddgsim),2),float)
    ddGsim[:,0] = ddgsim
    ddGsim[:,1] = ddgsim_err
    ddGexp = np.zeros((len(ddgexp),2),float)
    ddGexp[:,0] = ddgexp
    ddGexp[:,1] = ddgexp_err

    n_contacts = len(contacts)

    return epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, whamFree, state_labels, state_bounds, ddGsim, ddGexp

def plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=False):
    """ Get  """
    maxeps = max(epsilons)
    mineps = min(epsilons)
    if iteration == 0:
        plt.pcolor(epsilon_map,cmap=plt.get_cmap("Blues"))
        cbar = plt.colorbar()
        for i in range(len(contacts)):
            plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map(1.0))
    else:
        plt.pcolor(epsilon_map,cmap=plt.get_cmap("Blues"))
        cbar = plt.colorbar()
        cbar.set_clim(mineps,maxeps)
        for i in range(len(contacts)):
            plt.plot(contacts[i,0]-0.5,contacts[i,1]-0.5,marker='s',ms=6,markeredgecolor="k",markerfacecolor=new_map((epsilons[i]-mineps)/maxeps))

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
    plt.xlabel("$\\Delta\\Delta G_{sim}$")
    plt.ylabel("$\\Delta\\Delta G_{exp}$")
    if individual:
        plt.title("%s iteration %d" % (name,iteration))
    plt.xlim(-.2,7)
    plt.ylim(-.2,7)
    plt.grid(True)
    lg = plt.legend(loc=4)
    lg.draw_frame(False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    epsilons, epsilon_map, n_residues, contacts, n_contacts, Tf, whamFree, state_labels, state_bounds, ddGsim, ddGexp = get_iteration_data(name,iteration)

    plt.figure()
    plot_free_energy(name,iteration,n_contacts,Tf,whamFree,state_labels,state_bounds,individual=True)
    plt.savefig("%s/Mut_%d/FreeEnergy_Q.pdf" % (name,iteration))
    plt.close()

    plt.figure()
    plot_epsilon_map(name,iteration,epsilons,epsilon_map,contacts,n_residues,individual=True)
    plt.savefig("%s/Mut_%d/current_epsilon_map.pdf" % (name,iteration))
    plt.close()

    plt.figure()
    plot_epsilon_histogram(name,iteration,epsilons,individual=True)
    plt.savefig("%s/Mut_%d/current_epsilon_hist.pdf" % (name,iteration))
    plt.close()

    plt.figure()
    plot_ddG_comparison(name,iteration,ddGsim,ddGexp,individual=True)
    plt.savefig("%s/Mut_%d/compareddG.pdf" % (name,iteration))
    plt.close()

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
    plt.savefig("%s/Mut_%d/summary_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.show()
