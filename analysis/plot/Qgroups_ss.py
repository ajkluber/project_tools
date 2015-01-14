import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import argparse

from project_tools.analysis.plot.QivsQ import get_some_iteration_data
from project_tools.analysis.plot.drawss import add_secondary_struct_icons

def get_secondary_structure_Qgroups(name,contacts,n_contacts):
    """ Get Qgroups based on secondary structural elements """ 

    if not os.path.exists("%s/Qgroups_ss" % name):
        os.mkdir("%s/Qgroups_ss" % name)
    if os.path.exists("%s/Qgroups_ss/colors.txt" % name):
        colors = [ x.rstrip("\n") for x in open("%s/Qgroups_ss/colors.txt" % name,"r").readlines() ] 
    else: 
        colors = []

    ## Read in secondary structure assignment.
    ss_labels = []
    ss_bounds = []
    for line in open("%s/secondary_structure.txt" % name,"r").readlines():
        ss_labels.append(line.split()[0])
        ss_bounds.append([int(line.split()[1]),int(line.split()[2])])
    n_ss_elements = len(ss_labels)
    n_ss_groups = n_ss_elements + n_ss_elements*(n_ss_elements - 1)/2

    if not os.path.exists("%s/Qgroups_ss/group0.dat" % name):
        ## Group contacts as between secondary structural elements. 
        ## For each contact determine the group it belongs.
        Qgrp_conts = []
        Qgrp_indxs = []
        for n in range(n_ss_elements):
            temp = [ [] for m in range(n,n_ss_elements) ]
            temp2 = [ [] for m in range(n,n_ss_elements) ]
            Qgrp_conts.append(temp) 
            Qgrp_indxs.append(temp2)

        for i in range(n_contacts):
            cont = contacts[i]
            for n in range(n_ss_elements):
                if (cont[0] >= ss_bounds[n][0]) and (cont[0] <= ss_bounds[n][1]):
                    for m in range(n,n_ss_elements):
                        if (cont[1] >= ss_bounds[m][0]) and (cont[1] <= ss_bounds[m][1]):
                            Qgrp_conts[n][m-n].append(list(cont))
                            Qgrp_indxs[n][m-n].append(i)
                        else:
                            continue
                else:
                    continue

        Qgrp_indxs, colors = plot_Qgroups_ss_map(name,Qgrp_conts,Qgrp_indxs,n_ss_elements,ss_labels,ss_bounds,colors)
    else:
        n_groups = len(open("%s/Qgroups_ss/labels.txt" % name, "r").readlines())
        Qgrp_indxs = [ np.loadtxt("%s/Qgroups_ss/group%d.dat" % (name,x),dtype=int) for x in range(n_groups) ]
        colors = [ x.rstrip("\n") for x in open("%s/Qgroups_ss/colors.txt" % name, "r").readlines() ]
    
    return Qgrp_indxs, colors, ss_labels, ss_bounds

def plot_Qgroups_ss_map(name,Qgrp_conts,Qgrp_indxs,n_ss_elements,ss_labels,ss_bounds,colors):
    """ Plot the ss Qgroups on a contact map and save Qgroups"""

    if colors == []:
        needcolors = True
    else:
        needcolors = False

    colornames = matplotlib.colors.cnames.keys()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    i = 0
    group_name = ""
    colorstring = ""
    labelstring = ""
    temp_Qgroups = []
    color_indx = 0
    for n in range(n_ss_elements):
        for m in range(len(Qgrp_conts[n])):
            if not len(Qgrp_conts[n][m]) == 0:
                np.savetxt("%s/Qgroups_ss/group%d.dat" % (name,i), np.array(Qgrp_indxs[n][m]), fmt="%d")
                
                temp_Qgroups.append(np.array(Qgrp_indxs[n][m]).astype(int))
                labelstring += "%d %d\n" % (n,m)
                if needcolors:
                    somecolor = colornames[np.random.randint(len(colornames))]
                    while (somecolor in colors):
                        somecolor = colornames[np.random.randint(len(colornames))]
                    colorstring += somecolor + "\n"
                    colors.append(somecolor)
                else:
                    somecolor = colors[color_indx]
                color_indx += 1
                
                for p in range(len(Qgrp_conts[n][m])):
                    if p == 0:
                        print "group %d %d plotted in color %s" % (n,m,somecolor)
                        #ax.plot(Qgrp_conts[n][m][p][0],Qgrp_conts[n][m][p][1],marker='s',color=somecolor,label="%d %d" % (n,m))
                        ax.plot(Qgrp_conts[n][m][p][0],Qgrp_conts[n][m][p][1],marker='s',color=somecolor)
                    else:
                        ax.plot(Qgrp_conts[n][m][p][0],Qgrp_conts[n][m][p][1],marker='s',color=somecolor)
                i += 1
            else:
                pass

    if needcolors:
        open("%s/Qgroups_ss/colors.txt" % name, "w").write(colorstring)
        
    open("%s/Qgroups_ss/labels.txt" % name, "w").write(labelstring)
    ticks = []
    for a,b in ss_bounds:
        ticks.append(a)

    N = int(max(contacts.ravel()))
    add_secondary_struct_icons(ax,ss_labels,ss_bounds)
    plt.xlim(0,N)
    plt.ylim(0,N)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_yticklabels(), visible=False)
    plt.grid(True)
    #plt.legend(loc=4)
    plt.title("Qgroups for %s" % name)
    plt.savefig("%s/Qgroups_ss/group_map.png" % name)
    plt.savefig("%s/Qgroups_ss/group_map.pdf" % name)
    plt.show()
    return temp_Qgroups, colors
    
def plot_secondary_structure_Qgroups(name,iteration,Qbins,Qi_vs_Q,n_bins,contacts,n_contacts,state_labels,state_bounds,Qgrp_indxs,colors):

    plt.figure(1)
    plt.figure(2)
    n_grps = len(Qgrp_indxs)
    Qgrp_vs_Q = np.zeros((n_bins,n_grps),float)
    Qgrp_vs_Q_nrm = np.zeros((n_bins,n_grps),float)
    for n in range(n_grps):
        if Qgrp_indxs[n].shape == ():
            Qgrp_vs_Q[:,n] = Qi_vs_Q[:,Qgrp_indxs[n]]
            Qgrp_vs_Q_nrm[:,n] = Qi_vs_Q[:,Qgrp_indxs[n]]
        else:
            n_grp_members = Qgrp_indxs[n].shape[0]
            Qgrp_vs_Q[:,n] = sum(Qi_vs_Q[:,Qgrp_indxs[n]].T)
            Qgrp_vs_Q_nrm[:,n] = sum(Qi_vs_Q[:,Qgrp_indxs[n]].T)/n_grp_members

        plt.figure(1)
        plt.plot(Qbins,Qgrp_vs_Q[:,n],color=colors[n],lw=2)
        plt.figure(2)
        plt.plot(Qbins,Qgrp_vs_Q_nrm[:,n],color=colors[n],lw=2)

    plt.figure(1)
    plt.xlim(0,max(Qbins))
    #plt.ylim(0,1)
    plt.title("$Q_{ss}$ unnormalized %s iteration %d" % (name, iteration),fontsize=20)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{ss}$",fontsize=20)
    plt.savefig("%s/iteration_%d/plots/QssvsQ_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QssvsQ_%s_%d.png" % (name,iteration,name,iteration))

    plt.figure(2)
    plt.xlim(0,max(Qbins))
    plt.ylim(0,1)
    plt.title("$Q_{ss}$ %s iteration %d" % (name, iteration),fontsize=20)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{ss}$",fontsize=20)
    plt.savefig("%s/iteration_%d/plots/QssvsQ_nrm_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QssvsQ_nrm_%s_%d.png" % (name,iteration,name,iteration))
    plt.show()


if __name__ == '__main__':
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

    print "  getting Qgroups:"
    Qgrp_indxs, colors, ss_labels, ss_bounds = get_secondary_structure_Qgroups(name,contacts.astype(int),n_contacts)

    print "  plotting Qgroup_ss vs Q:"
    #print Qgrp_indxs, colors
    plot_secondary_structure_Qgroups(name,iteration,Qbins,Qi_vs_Q,n_bins,contacts,n_contacts,state_labels,state_bounds,Qgrp_indxs,colors)
