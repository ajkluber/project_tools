import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import argparse

from project_tools.analysis.plot.plot_Qgroups import get_some_iteration_data


def plot_secondary_structure_Qgroups(name,iteration,Qbins,Qi_vs_Q,n_bins,contacts,n_contacts,state_labels,state_bounds):

    ## Read in secondary structure assignment.
    ss_labels = []
    ss_bounds = []
    for line in open("%s/secondary_structure.txt" % name,"r").readlines():
        ss_labels.append(line.split()[0])
        ss_bounds.append([int(line.split()[1]),int(line.split()[2])])

    ## For N ss elements, there are N + N(N-1)/2 groups.
    n_ss_elements = len(ss_labels)
    n_ss_groups = n_ss_elements + (n_ss_elements*(n_ss_elements - 1)/2)
    #print n_ss_elements          ## DEUBUGGING
    #print n_ss_groups          ## DEUBUGGING
    
    Qgroups = []
    for n in range(n_ss_elements):
        temp = [ [] for m in range(n,n_ss_elements) ]
        #print temp             ## DEBUGGING
        Qgroups.append(temp) 
    
    unsorted = []
    #print Qgroups          ## DEUBUGGING
    

    ## Group contacts as between secondary structural elements. 
    ## For each contact determine the group it belongs.
    checker = []
    for i in range(n_contacts):
        cont = contacts[i]
        flag1 = [ (cont[0] >= ss_bounds[x][0]) and (cont[0] < ss_bounds[x][1]) for x in range(n_ss_elements) ]
        flag2 = [ (cont[1] >= ss_bounds[x][0]) and (cont[1] < ss_bounds[x][1]) for x in range(n_ss_elements) ]
        if (not any(flag1)) or (not any(flag2)): 
            unsorted.append(cont)
        else:
            grp1 = flag1.index(True)
            grp2 = flag2.index(True) - grp1
            Qgroups[grp1][grp2].append(cont)

    print Qgroups
    print unsorted

    raise SystemExit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    parser.add_argument('--n_bins', type=int, default=35, help='Number of bins along Q.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    n_bins = args.n_bins

    ## Get some iteration data
    epsilons, loops, n_residues, contacts, n_contacts, Tf, state_labels, state_bounds, Qbins, Qi_vs_Q = get_some_iteration_data(name,iteration,n_bins)


    print "  plotting Qgroups:"
    plot_secondary_structure_Qgroups(name,iteration,Qbins,Qi_vs_Q,n_bins,contacts,n_contacts,state_labels,state_bounds)
