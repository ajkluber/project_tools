""" Library for plotting free energy surfaces (or pmfs). 

Description:

    Functions that read output data files from WHAM calculations and
plot 1D or 2D free energy surfaces (aka potentials of mean force (1)). 


To Do:

- Create command line functionality that can be called from a PBS script.
- Write a function that submits PBS scripts to plot quantities.


References:
(1) http://en.wikipedia.org/wiki/Potential_of_mean_force
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
import argparse
import subprocess as sb

def get_data(coord):
    if coord in ["Q","Qh","Qnh"]:
        #x = np.loadtxt(coord+"prob.dat")
        x = np.loadtxt(coord+".dat")
        x /= max(x)
    elif coord == "Nh":
        x = np.loadtxt(coord+".dat") 
        x /= max(x)
    elif coord == "rmsd":
        dummy, x = np.loadtxt(coord+".xvg",unpack=True)
    elif coord == "Rg":
        dummy, x = np.loadtxt("radius_cropped.xvg",unpack=True)
    else:
        print "ERROR!"
        print "  Coordinate: ",coord," not found"
    return x

def plot_aggregated_data(System,append_log):
    ''' Plot 1D and 2D pmfs for aggregated equilibrium simulations.'''
    
    append_log(System.subdir,"Starting: Plotting_Agg_Data")
    cwd = os.getcwd()
    sub = System.subdir+"/"+System.mutation_active_directory
    os.chdir(cwd+"/"+sub)
    temps = [ x.split('_')[0] for x in open("T_array.txt","r").readlines() ] 
    unique_temps = []
    counts = []
    for t in temps:
        if t not in unique_temps:
            unique_temps.append(t)
            counts.append(temps.count(t))
        else:
            pass

    cwd2 = os.getcwd()
    coords = ["Q","Qh","Qnh","Nh","Rg","rmsd"]
    coord_pairs = [("Qh","Q"),("Qh","Qnh"),("Qnh","Q"),("Nh","Qnh"),("Rg","rmsd")]

    for i in range(len(unique_temps)):
        
        T = unique_temps[i]
        os.chdir(T+"_agg")
        print "  Plotting PMFs for ",T

        for crd in coords:
            if not os.path.exists("pmfs/"+crd+"_pmf.pdf"):
                try:
                    plot_1D_pmf(crd,System.subdir+" "+T)
                    print "    Plotted pmf ",crd
                except:
                    print "    Plotting pmf ",crd, " ** didn't work **. skipping."

        for crd1,crd2 in coord_pairs:
            if not os.path.exists("pmfs/"+crd1+"_"+crd2+"_pmf.pdf"):
                try:
                    plot_2D_pmf(crd1,crd2,System.subdir+" "+T)
                    print "    Plotted 2D pmf ",crd1, " vs ", crd2
                except:
                    print "    Plotting 2D pmf ",crd1, " vs ", crd2," ** didn't work **. skipping."

        os.chdir(cwd2)

    os.chdir(cwd)

    append_log(System.subdir,"Finished: Plotting_Agg_Data")

def plot_1D_pmf(coord,title):
    ''' Plot a 1D pmf for a coordinate.'''

    x = get_data(coord)

    path = os.getcwd()
    savedir = path+"/pmfs"
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    if coord in ["Rg","rmsd"]:
        skip = 80
    else:
        skip = 4 
    vals = np.unique(list(x))[::skip]

    n,bins = np.histogram(x,bins=vals,density=True)
    np.savetxt(savedir+"/"+coord+"_n.dat",n,delimiter=" ",fmt="%.4f")
    np.savetxt(savedir+"/"+coord+"_bins.dat",bins,delimiter=" ",fmt="%.4f")

    pmf = -np.log(n)
    pmf -= min(pmf)

    plt.figure()
    plt.plot(bins[1:]/max(bins),pmf)
    plt.xlabel(coord,fontsize="xx-large")
    plt.ylabel("F("+coord+") / kT",fontsize="xx-large")
    plt.title("F("+coord+") "+title,fontsize="xx-large")
    plt.ylim(0,6)
    plt.xlim(0,1)
    plt.savefig(savedir+"/"+coord+"_pmf.pdf")

def plot_2D_pmf(coord1,coord2,title=""):

    path = os.getcwd()
    savedir = path+"/pmfs"
    if os.path.exists(savedir) == False:
        os.mkdir(savedir)

    x = get_data(coord1)
    y = get_data(coord2)

    if (coord1 in ["Rg","rmsd"]) or (coord2 in ["Rg","rmsd"]):
        skip = 80
    else:
        skip = 4 
    xvals = np.unique(list(x))[::skip]
    yvals = np.unique(list(y))[::skip]

    hist,xedges,yedges = np.histogram2d(x,y,bins=[xvals,yvals],normed=True)
    hist[hist == 0.0] == 1.0
    pmf = -np.log(hist)
    pmf -= min(pmf.ravel())
    X,Y = np.meshgrid(yedges[:-1],xedges[:-1])

    levels = np.arange(0,6,0.5)

    plt.figure()
    plt.subplot(1,1,1,aspect=1)
    plt.contourf(X,Y,pmf,levels=levels)
    if (coord1.startswith("Q")) or (coord2.startswith("Q")) or (coord1 == "Nh") or (coord2 == "Nh"):
        plt.xlim(0,1)
        plt.ylim(0,1)
    else: 
        plt.xlim(0,max(x))
        plt.ylim(0,max(y))
    plt.xlabel(coord2,fontsize="xx-large")
    plt.ylabel(coord1,fontsize="xx-large")
    plt.title("F("+coord1+","+coord2+") "+title,fontsize="xx-large")
    cbar = plt.colorbar()
    cbar.set_label("F / kT",fontsize="xx-large")
    #plt.savefig(path+"/Q_Qh_"+T+"_pmf.pdf")
    plt.savefig(savedir+"/"+coord1+"_"+coord2+"_pmf.pdf")
