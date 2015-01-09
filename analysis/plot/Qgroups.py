import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import argparse

from project_tools.analysis.plot.QivsQ import get_some_iteration_data

def get_groupset_groups(name,groupset):
    """ Get a user defined set of Qgroups from file """

    colornames = matplotlib.colors.cnames.keys()
    if not os.path.exists("%s/Qgroups_%d" % (name,groupset)):
        print "ERROR! %s/Qgroups_%d DOES NOT EXIST!" % (name,groupset)
        print " Specify sets of groups for groupset %d " % groupset
        print "Exiting"
        raise SystemExit


    Qgrp_indxs = []
    grp_labels = [ x.rstrip("\n") for x in open("%s/Qgroups_%d/labels.txt" % (name,groupset),"r").readlines() ] 
    n_grps = len(grp_labels)
    for i in range(n_grps):
        Qgrp_indxs.append(np.loadtxt("%s/Qgroups_%d/group%d.dat" % (name,groupset,i),dtype=int))

    if not os.path.exists("%s/Qgroups_%d/colors.txt" % (name,groupset)):
        colors = []
        colorstring = ""
        for k in range(n_grps):
            somecolor = colornames[np.random.randint(len(colornames))] 
            colorstring += somecolor + "\n"
            colors.append(somecolor)
        open("%s/Qgroups_%d/colors.txt" % (name,groupset), "w").write(colorstring)
    else:
        colors = [ x.rstrip("\n") for x in open("%s/Qgroups_%d/colors.txt" % (name,groupset),"r").readlines() ] 

    ## TO DO:
    ## Plot Group map legend in <name>/Qgroups_<#>/
    if not os.path.exists("%s/Qgroups_%d/group_map.png" % (name,groupset)):
        contacts = np.loadtxt("%s/contacts.dat" % name,dtype=int)
        plt.figure()
        for n in range(n_grps):
            qgrp = Qgrp_indxs[n]
            for i in range(len(qgrp)):
                indx = qgrp[i]
                if i == 0:
                    plt.plot(contacts[indx,0],contacts[indx,1],marker='s',color=colors[n],label=grp_labels[n])
                else:
                    plt.plot(contacts[indx,0],contacts[indx,1],marker='s',color=colors[n])
        plt.legend(loc=4)
        N = max(contacts.ravel())
        plt.xticks(range(0,N,10))
        plt.yticks(range(0,N,10))
        plt.grid(True)
        plt.title("$Q_{groups}$ groupset %d for %s" % (groupset, name))
        plt.savefig("%s/Qgroups_%d/group_map.png" % (name,groupset))
        plt.savefig("%s/Qgroups_%d/group_map.pdf" % (name,groupset))

    return Qgrp_indxs, n_grps, grp_labels, colors

def plot_Qgroups_vs_Q(name,iteration,groupset,Qgrp_indxs,n_grps,grp_labels,colors,Qbins,Qi_vs_Q):
    """ Plot Qgroups vs Q for a groupset """

    plt.figure(1)
    plt.figure(2)
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
        plt.plot(Qbins,Qgrp_vs_Q[:,n],color=colors[n],lw=2,label=grp_labels[n])
        plt.figure(2)
        plt.plot(Qbins,Qgrp_vs_Q_nrm[:,n],color=colors[n],lw=2,label=grp_labels[n])

    plt.figure(1)
    plt.xlim(0,max(Qbins))
    #plt.ylim(0,1)
    plt.title("$Q_{group}$ groupset %d %s iteration %d" % (groupset,name,iteration),fontsize=18)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{group}$",fontsize=20)
    lg = plt.legend(loc=2)
    lg.draw_frame(False)
    plt.savefig("%s/iteration_%d/Qgroups_%d/QgroupvsQ_%s_%d.pdf" % (name,iteration,groupset,name,iteration))
    plt.savefig("%s/iteration_%d/Qgroups_%d/QgroupvsQ_%s_%d.png" % (name,iteration,groupset,name,iteration))

    plt.figure(2)
    plt.xlim(0,max(Qbins))
    plt.ylim(0,1)
    plt.title("$Q_{group}$ groupset %d %s iteration %d" % (groupset,name,iteration),fontsize=18)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{group}$",fontsize=20)
    lg = plt.legend(loc=2)
    lg.draw_frame(False)
    plt.savefig("%s/iteration_%d/Qgroups_%d/QgroupvsQ_nrm_%s_%d.pdf" % (name,iteration,groupset,name,iteration))
    plt.savefig("%s/iteration_%d/Qgroups_%d/QgroupvsQ_nrm_%s_%d.png" % (name,iteration,groupset,name,iteration))
    plt.show()

    return Qgrp_vs_Q,Qgrp_vs_Q_nrm

def multiplot_Qgroups_vs_Q(name,iterations,n_contacts,groupset,n_grps,grp_labels,colors,Qbins,All_QgrpsvsQ,All_QgrpsvsQ_nrm):
    """ Plot a set of iterations all together on one multiplot """

    n_iters = len(iterations)

    maxy = max([ max(All_QgrpsvsQ[i][-1,:]) for i in range(n_iters) ]) 
    savestring = ""

    fig, axes = plt.subplots(1,n_iters,sharey=True,figsize=(5*n_iters,5.5))

    for m in range(n_iters):
        savestring += "_%d" % iterations[m]
        ax = axes[m]
        ax.set_title("iteration %d" % iterations[m])
        for n in range(n_grps):
            ax.plot(Qbins,All_QgrpsvsQ[m][:,n],color=colors[n],lw=2,label=grp_labels[n])
        lg = ax.legend(loc=2)
        lg.draw_frame(False)
        
    
    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("$Q_{group}$",fontsize=20)
    fig.text(0.5,0.04,"Foldedness $Q$",fontsize=20,ha='center',va='center')
    fig.text(0.5,0.96,"$Q_{groups}$ vs $Q$ for groupset %d %s" % (groupset,name),fontsize=18,ha='center',va='center')
    fig.savefig("%s/Qgroups_%d/QgrpvsQ%s.pdf" % (name,groupset,savestring))
    fig.savefig("%s/Qgroups_%d/QgrpvsQ%s.png" % (name,groupset,savestring))
    plt.clf()


    fig, axes = plt.subplots(1,n_iters,sharey=True,figsize=(5*n_iters,5.5))
    for p in range(n_iters):
        ax = axes[p]
        ax.set_title("iteration %d" % iterations[p])

        for n in range(n_grps):
            ax.plot(Qbins,All_QgrpsvsQ_nrm[p][:,n],color=colors[n],lw=2,label=grp_labels[n])
        lg = ax.legend(loc=2)
        lg.draw_frame(False)
    
    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("$Q_{group}$",fontsize=20)
    fig.text(0.5,0.04,"Foldedness $Q$",fontsize=20,ha='center',va='center')
    fig.text(0.5,0.96,"$Q_{groups}$ vs $Q$ for groupset %d %s" % (groupset,name),fontsize=18,ha='center',va='center')
    fig.savefig("%s/Qgroups_%d/QgrpvsQ_nrm%s.pdf" % (name,groupset,savestring))
    fig.savefig("%s/Qgroups_%d/QgrpvsQ_nrm%s.png" % (name,groupset,savestring))
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iterations', type=int, nargs="+", required=True, help='Iteration to plot.')
    parser.add_argument('--groupset', type=int, required=True, help='Set of groups to plot.')
    parser.add_argument('--n_bins', type=int, default=35, help='Number of bins along Q.')
    args = parser.parse_args()

    name = args.name
    groupset = args.groupset
    iterations = args.iterations
    n_bins = args.n_bins
    if len(iterations) == 1:
        multiplot = False
    else:
        multiplot = True

    ## Get  group set. Qgroup
    print "  getting Qgroups for groupset %d" % groupset
    Qgrp_indxs, n_grps, grp_labels, colors = get_groupset_groups(name,groupset)

    
    raise SystemExit

    All_QgrpsvsQ = []
    All_QgrpsvsQ_nrm = []
    
    savestring = ""
    for n in range(len(iterations)):
        iteration = iterations[n]
        savestring += "_%d" % iteration
        if not os.path.exists("%s/iteration_%d/Qgroups_%d" % (name,iteration,groupset)):
            os.mkdir("%s/iteration_%d/Qgroups_%d" % (name,iteration,groupset))

        ## Get some iteration data
        epsilons, loops, n_residues, contacts, n_contacts, Tf, state_labels, state_bounds, Qbins, Qi_vs_Q = get_some_iteration_data(name,iteration,n_bins)

        ## Plot Qgroups vs Q
        print "  plotting Qgroups vs Q for groupset %d" % groupset
        print "    saving as: %s/iteration_%d/Qgroups_%d/QgroupvsQ_nrm_%s_%d.pdf" % (name,iteration,groupset,name,iteration)
        Qgrp_vs_Q, Qgrp_vs_Q_nrm = plot_Qgroups_vs_Q(name,iteration,groupset,Qgrp_indxs,n_grps,grp_labels,colors,Qbins,Qi_vs_Q)
        All_QgrpsvsQ.append(Qgrp_vs_Q)
        All_QgrpsvsQ_nrm.append(Qgrp_vs_Q_nrm)

    if multiplot:
        print "  plotting multiplot for all iterations"
        print "    saving as: %s/Qgroups_%d/QgrpvsQ_nrm%s.png(pdf)" % (name,groupset,savestring)
        multiplot_Qgroups_vs_Q(name,iterations,n_contacts,groupset,n_grps,grp_labels,colors,Qbins,All_QgrpsvsQ,All_QgrpsvsQ_nrm)
