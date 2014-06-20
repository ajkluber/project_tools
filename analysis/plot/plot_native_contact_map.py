#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import argparse

import mdtraj as md

def plot_native_state_contact_map(title):

    colors = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

    if os.path.exists("contact_pairs.dat") and os.path.exists("contact_probabilities.dat"):
        pairs = np.loadtxt("contact_pairs.dat")
        probability = np.loadtxt("contact_probabilities.dat")
    else:
        print "  Loading BeadBead.dat"
        beadbead = np.loadtxt("BeadBead.dat",dtype=str)
        sigij = beadbead[:,5].astype(float)
        epsij = beadbead[:,6].astype(float)
        deltaij = beadbead[:,7].astype(float)
        interaction_numbers = beadbead[:,4].astype(str)
        pairs = beadbead[:,:2].astype(int)
        pairs -= np.ones(pairs.shape,int)
        np.savetxt("contact_pairs.dat",pairs)

        print "  Computing distances with mdtraj..."
        traj = md.load("traj.xtc",top="Native.pdb")
        distances = md.compute_contacts(traj,pairs)
        contacts = (distances[0][:] <= 1.2*sigij).astype(int)
        print "  Computing contact probability..."
        probability = sum(contacts.astype(float))/contacts.shape[0]
        np.savetxt("contact_probabilities.dat",probability)

    Qref = np.loadtxt("Qref_cryst.dat")
    C = np.zeros(Qref.shape,float)

    for k in range(len(pairs)):
        C[pairs[k][0],pairs[k][1]] = probability[k]

    print "  Plotting..."
    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    ax = plt.subplot(1,1,1,aspect=1)
    plt.pcolor(C,cmap=new_map)
    for k in range(len(pairs)):
        if probability[k] > 0.01:
            plt.plot(pairs[k][1],pairs[k][0],marker='s',ms=3.0,markeredgecolor=new_map(probability[k]),color=new_map(probability[k]))
        else:
            continue
    plt.xlim(0,len(Qref))
    plt.ylim(0,len(Qref))
    #plt.text(10,70,name.upper(),fontsize=70,color="r")
    ax = plt.gca()
    cbar = plt.colorbar()
    cbar.set_clim(0,1)
    cbar.set_label("Contact probability",fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    plt.xlabel("Residue i",fontsize=20)
    plt.ylabel("Residue j",fontsize=20)
    #plt.title("Native State Contact Map "+title,fontsize=20)
    plt.title(title)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(15)
    print "  Saving..."
    plt.savefig("native_state_contact_map.pdf")
    #plt.show()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Perform WHAM.')
    parser.add_argument('--title', type=str, help='Method',required=True)
    args = parser.parse_args()
    print args

    temps = [ t[:-1] for t in open("T_array.txt","r").readlines() ]

    cwd = os.getcwd()
    for tdir in temps:
        print "Entering directory ",tdir
        os.chdir(tdir)
        if os.path.exists("native_state_contact_map.pdf"):
            print "  Skipping"
        else:
            title = "%s $T = %sK$" % (args.title,tdir.split("_")[0])
            plot_native_state_contact_map(title)
        os.chdir(cwd)
