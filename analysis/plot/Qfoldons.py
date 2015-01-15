import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import argparse

from project_tools.analysis.plot.QivsQ import get_contact_probability_versus_Q

def get_foldon_groups(name):
    """ Get a user defined set of Qfoldons from file """
    if not os.path.exists("%s/Qfoldons" % name):
        print "ERROR! %s/Qfoldons DOES NOT EXIST!" % name
        print " Specify sets of contacts"
        print "Exiting"
        raise SystemExit

    labels = [ x.rstrip("\n") for x in open("%s/Qfoldons/labels.txt" % name,"r").readlines() ] 
    colors = [ x.rstrip("\n") for x in open("%s/Qfoldons/colors.txt" % name,"r").readlines() ] 
    n_grps = len(colors)

    Qgrp_indxs = []
    for i in range(n_grps):
        Qgrp_indxs.append(np.loadtxt("%s/Qfoldons/group%d.dat" % (name,i),dtype=int))

    return Qgrp_indxs, n_grps, colors, labels

def plot_Qfoldons_vs_Q(name,iteration,Qgrp_indxs,n_grps,colors,Qbins,Qi_vs_Q,labels):
    """ Plot Qfoldons vs Q"""

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
        plt.plot(Qbins,Qgrp_vs_Q[:,n],color=colors[n],lw=2,label=labels[n])
        plt.figure(2)
        plt.plot(Qbins,Qgrp_vs_Q_nrm[:,n],color=colors[n],lw=2,label=labels[n])

    plt.figure(1)
    plt.xlim(0,max(Qbins))
    #plt.ylim(0,1)
    plt.title("$Q_{group}$ %s iteration %d" % (name,iteration),fontsize=18)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{group}$",fontsize=20)
    lg = plt.legend(loc=2)
    lg.draw_frame(False)
    plt.savefig("%s/iteration_%d/plots/QfoldonsvsQ_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QfoldonsvsQ_%s_%d.png" % (name,iteration,name,iteration))

    plt.figure(2)
    plt.xlim(0,max(Qbins))
    plt.ylim(0,1)
    plt.title("$Q_{group}$ %s iteration %d" % (name,iteration),fontsize=18)
    plt.xlabel("$Q$",fontsize=20)
    plt.ylabel("$Q_{group}$",fontsize=20)
    lg = plt.legend(loc=2)
    lg.draw_frame(False)
    plt.savefig("%s/iteration_%d/plots/QfoldonsvsQ_nrm_%s_%d.pdf" % (name,iteration,name,iteration))
    plt.savefig("%s/iteration_%d/plots/QfoldonsvsQ_nrm_%s_%d.png" % (name,iteration,name,iteration))
    plt.show()

    return Qgrp_vs_Q,Qgrp_vs_Q_nrm

def multiplot_Qfoldons_vs_Q(name,iterations,n_grps,colors,Qbins,All_QgrpsvsQ,All_QgrpsvsQ_nrm,labels):
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
            ax.plot(Qbins,All_QgrpsvsQ[m][:,n],color=colors[n],lw=2,label=labels[n])
        lg = ax.legend(loc=2)
        lg.draw_frame(False)
        
    
    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("$Q_{group}$",fontsize=20)
    fig.text(0.5,0.04,"Foldedness $Q$",fontsize=20,ha='center',va='center')
    fig.text(0.5,0.96,"$Q_{groups}$ vs $Q$ for %s" % name,fontsize=18,ha='center',va='center')
    fig.savefig("%s/plots/QfoldonsvsQ%s.pdf" % (name,savestring))
    fig.savefig("%s/plots/QfoldonsvsQ%s.png" % (name,savestring))
    plt.clf()


    fig, axes = plt.subplots(1,n_iters,sharey=True,figsize=(5*n_iters,5.5))
    for p in range(n_iters):
        ax = axes[p]
        ax.set_title("iteration %d" % iterations[p])

        for n in range(n_grps):
            ax.plot(Qbins,All_QgrpsvsQ_nrm[p][:,n],color=colors[n],lw=2,label=labels[n])
        lg = ax.legend(loc=2)
        lg.draw_frame(False)
    
    plt.subplots_adjust(wspace=0)
    axes[0].set_ylabel("$Q_{group}$",fontsize=20)
    fig.text(0.5,0.04,"Foldedness $Q$",fontsize=20,ha='center',va='center')
    fig.text(0.5,0.96,"$Q_{groups}$ vs $Q$ for %s" % name,fontsize=18,ha='center',va='center')
    fig.savefig("%s/plots/QfoldonsvsQ_nrm%s.pdf" % (name,savestring))
    fig.savefig("%s/plots/QfoldonsvsQ_nrm%s.png" % (name,savestring))
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iterations', type=int, nargs='+',required=True, help='Iteration to plot.')
    parser.add_argument('--n_bins', type=int, default=35, help='Number of bins along Q.')
    args = parser.parse_args()

    name = args.name
    iterations = args.iterations
    n_bins = args.n_bins
    if len(iterations) == 1:
        multiplot = False
    else:
        multiplot = True


    ## Get  group set. 
    print "  getting Qfoldons for "
    Qgrp_indxs, n_grps, colors, labels = get_foldon_groups(name)
    
    All_QgrpsvsQ = []
    All_QgrpsvsQ_nrm = []
    savestring = ""
    for n in range(len(iterations)):
        iteration = iterations[n]
        savestring += "_%d" % iteration
        if not os.path.exists("%s/iteration_%d/plots" % (name,iteration)):
            os.mkdir("%s/iteration_%d/plots" % (name,iteration))

        ## Get Qi_vsQ
        Qbins, Qi_vs_Q = get_contact_probability_versus_Q(name,iteration,n_bins)

        ## Plot Qfoldons vs Q
        print "  plotting Qfoldons vs Q"
        print "    saving as: %s/iteration_%d/plots/QfoldonsvsQ_nrm_%s_%d.pdf" % (name,iteration,name,iteration)
        Qgrp_vs_Q, Qgrp_vs_Q_nrm = plot_Qfoldons_vs_Q(name,iteration,Qgrp_indxs,n_grps,colors,Qbins,Qi_vs_Q,labels)
        All_QgrpsvsQ.append(Qgrp_vs_Q)
        All_QgrpsvsQ_nrm.append(Qgrp_vs_Q_nrm)

    if multiplot:
        print "  plotting multiplot for all iterations"
        print "    saving as: %s/plots/QgrpvsQ_nrm%s.png(pdf)" % (name,savestring)
        multiplot_Qfoldons_vs_Q(name,iterations,n_grps,colors,Qbins,All_QgrpsvsQ,All_QgrpsvsQ_nrm,labels)
