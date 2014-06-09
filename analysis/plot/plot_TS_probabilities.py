""" Plot the probabilities of forming native contacts in the transition state.


Description:

    The transition state ensemble for a two-state folding protein is the 
ensemble of structures found at the top of the free energy barrier (when
an appropriate reaction coordinate is chosen).

"""

import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import mdtraj as md

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--Qref', type=str, help='Contact Map',required=True)
    parser.add_argument('--coord', type=str, help='Reaction coordinate',required=True)
    parser.add_argument('--bounds', type=float, nargs='+', help='TS bounds along reaction coordinate',required=True)
    args = parser.parse_args()

    Qref = np.loadtxt(args.Qref)
    x = np.loadtxt(args.coord)
    bounds = args.bounds
    N = len(Qref)

    colors = [('white')] + [(cm.jet(i)) for i in xrange(1,256)]
    new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=256)

    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt("BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(int)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    pairs = pairs[ interaction_numbers != 0 ]
    sigij = sigij[ interaction_numbers != 0 ]

    print "  Computing distances with mdtraj..."
    traj = md.load("traj.xtc",top="Native.pdb")
    distances = md.compute_contacts(traj,pairs)
    contacts = (distances[0][:] <= 1.2*sigij).astype(int)
    
    keep_frames = (((x > bounds[0]).astype(int)*(x < bounds[1]).astype(int)) == 1)
    contacts = contacts[keep_frames,:]
     
    print "  Computing contact probability..."
    probability = sum(contacts.astype(float))/contacts.shape[0]

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
            plt.plot(pairs[k][1]+0.5,pairs[k][0]+0.5,marker='s',ms=3.0,markeredgecolor=new_map(probability[k]),color=new_map(probability[k]))
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
    plt.title("TS contact probability")
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(15)
    if not os.path.exists("plots"):
        os.makedirs("plots")
    print "  Saving..."
    plt.savefig("plots/TS_contact_probability.pdf")
    plt.show()
