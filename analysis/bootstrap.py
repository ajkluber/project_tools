""" Bootstrap error analysis for WHAM free energy curves.


Description:

    Bootstrap analysis can be done in any situation where no closed form
equations for the statistical error exist. Formulas for the error in WHAM
existence, but they require making assumptions about the data that may not be
true. Bootstrapping on the other hand is a general way to estimate the effects
of finite sampling.

    Bootstrapping involves using the observed data from simulations to concoct
many hypothetical observations. These concocted datasets are run through WHAM
then the variance of resulting outputs is taken as an estimate of the variance
of the final result from the observed data.


References:
(1) http://en.wikipedia.org/wiki/Bootstrapping_%28statistics%29
(2) http://www.alchemistry.org/wiki/Analyzing_Simulation_Results#Bootstrap_Sampling
(3) Efron, B.; Tibshirani, R. "An Introduction to the Bootstrap" (1993)
"""

import numpy as np
import subprocess as sb
import os
import time
import shutil
import argparse

import project_tools.analysis.wham as wham
import resample_histo     ## Cython extension 

def get_state_bounds(name,iteration):
    state_labels = []
    state_bounds = []
    for line in open("%s/iteration_%d/state_bounds.txt" % (name,iteration),"r"):
        state_labels.append(line.split()[0])
        state_bounds.append([int(line.split()[1]),int(line.split()[2])])
    return state_labels, state_bounds

def get_iteration_data():
    """ Get options from WHAM.config file """
    histos = []
    temps = []
    numBins = []
    start = []
    step = []

    for line in open("free.config","r"):
        if line.startswith("name"):
            h_name = line.split()[1]
            T = float(h_name.split("_")[1])
            temps.append(T)
            histos.append(np.loadtxt(h_name,dtype=float))

        if line.startswith("numBins"):
            numBins.append(int(line.split()[1]))
        if line.startswith("start"):
            start.append(float(line.split()[1]))
        if line.startswith("step"):
            step.append(float(line.split()[1]))

    n_samples = histos[0].shape[0]
    n_temps = len(histos)
    return histos, temps, numBins, start, step, n_temps, n_samples

## To perform bootstrapping need to run WHAM many times on datasets generated
## from simulation dataset. To be quick and memory efficient:
## Read in WHAM options from free.config -> Bins in energy and coordinate
## Read in Tf used for final output F(Q)
## Read in final output F(Q) curve 
## 
## Create bootstrap subdirectory, copy WHAM.jar executable and free.config (only output temperature Tf)
## into it and enter. 
## For n in N bootstrapping iterations:
##   1. Generate hypothetical datasets, save datasets.
##   2. Run WHAM
##   3. Copy output to filename unique to iteration
##   4. Delete hypothetical datasets (or let them be overwritten?)
##
## After WHAM iterations are done. Load in F(Q) from bootstrapping iterations
## Take the variance for each bin, save variance and plot.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    parser.add_argument('--num_boots', type=int, default=50, help='Number of bootstrapping iterations.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    n_boots = args.num_boots

    ## Temp used for parameter fitting
    Tuse = float(open("%s/iteration_%d/long_temps_last" % (name,iteration),"r").readlines()[0].split("_")[0])
    ## Folding temp
    #Tf = open("%s/Mut_%d/Tf.txt" % (name,iteration),"r").read().rstrip("\n")                       
    Tf = open("%s/iteration_%d/long_Tf" % (name,iteration),"r").read().rstrip("\n")                       
    startTF = Tuse
    deltaTF = 0.1
    ntempsF = 1
    state_labels, state_bounds = get_state_bounds(name,iteration)

    cwd = os.getcwd()
    #os.chdir("%s/Mut_%d/whamQ" % (name,iteration))
    os.chdir("%s/iteration_%d/long_wham" % (name,iteration))

    print "Getting free.config options..."
    histos, temps, numBins, start, step, n_temps, n_samples = get_iteration_data()

    n_Ebins = numBins[0]
    Ebins = np.arange(start[0],start[0] + (n_Ebins + 1)*step[0],step[0])
    n_Xbins = numBins[1]
    Xbins = np.arange(start[1],start[1] + (n_Xbins + 1)*step[1],step[1])

    if not os.path.exists("bootstrap"):
        os.mkdir("bootstrap")
    os.chdir("bootstrap")

    PROJECTS = os.environ["PROJECTS"]
    shutil.copy(PROJECTS+"/project_tools/analysis/WHAM.1.06.jar",".")

    ## Create WHAM input file boot.config
    print "Creating WHAM config file..."
    wham_config = wham.get_wham_config_basic(numBins[0],start[0],step[0],numBins[1],start[1],step[1],2)
    wham_config += "### list of histogram filenames and their temperatures ###\n"
    wham_config += "numFiles %d\n" % n_temps
    for i in range(n_temps):
        wham_config += "name hist_%d temp %.2f\n" % (i,temps[i])
    wham_config += "\n"
    wham_config += wham.get_wham_config_free_energy(round(startTF,1),deltaTF,ntempsF)
    open("free_boot.config","w").write(wham_config)

    print "Starting bootstrapping..."
    ## Run bootstrap iterations
    free_est = np.zeros((n_boots,numBins[1]-1),float)
    for n in range(n_boots):

        if os.path.exists("free_%d" % n):
            free_est[n,:] = np.loadtxt("free_%d" % n,usecols=(1,))
            print "  skip iteration %d" % n
            continue
        else:
            start_iter = time.clock()
            #def bootstrapping_iteration():
            for i in range(n_temps):
                ## Calculate the cumulative distribution for Energy.
                E = histos[i][:,0]
                Ebin_counts,Ebin_edges = np.histogram(E,bins=Ebins)
                Ecumul = np.array([ sum((Ebin_counts.astype(float)/float(n_samples))[:x]) for x in range(n_Ebins) ] + [1.0]).astype(np.double)
                Ebin_edges = Ebin_edges.astype(np.double)

                ## Calculate the cumulative distribution for some reaction coordinate.
                X = histos[i][:,1]
                Xbin_counts,Xbin_edges = np.histogram(X,bins=Xbins)
                Xcumul = np.array([ sum((Xbin_counts.astype(float)/float(n_samples))[:x]) for x in range(n_Xbins) ] + [1.0]).astype(np.double)
                Xbin_edges = Xbin_edges.astype(np.double)

                ## Generate sample data using cython function.
                #print "Generating histograms for temp %6.2f..." % temps[i]
                #starttime = time.clock()
                gen_data = resample_histo.generate_histo(n_samples,n_Ebins,Ebin_edges,Ecumul,n_Xbins,Xbin_edges,Xcumul)
                np.savetxt("hist_%d" % i,gen_data)      ## 
                #calctime = time.clock() - starttime
                #print "     Gen histo took: %f sec" % calctime

                ## DEBUGGING
                #gen_histo = np.zeros((Ebins),int)
                #gen_data = np.zeros((n_samples),float)
                #for i in range(n_samples):
                #    x = np.random.rand() 
                #    for n in range(1,Ebins):
                #        if (x > cumul[n-1]) and (x <= cumul[n]):
                #            gen_histo[n-1] = gen_histo[n-1] + 1
                #            gen_data[i] = 0.5*(bin_edges[n-1] + bin_edges[n])
                #            break
                #stoptime = time.clock()
                #print "Generating histogram took: %f sec" % (stoptime - starttime)
                #gen_histo,gen_bin_edges = np.histogram(gen_data,bins=bins)
                #plt.plot(0.5*(bin_edges[1:] + bin_edges[:-1]),bin_counts,lw=2,color='r')
                #plt.plot(0.5*(bin_edges[1:] + bin_edges[:-1]),gen_histo,lw=2,color='b')
                #plt.show()

            cmd1 = "java -jar WHAM.1.06.jar --config free_boot.config"      ## 'module load jdk' to use java
            #starttime = time.clock()
            sb.call(cmd1.split(),stdout=open("boot.out","w"),stderr=open("boot.err","w"))
            #calctime = time.clock() - starttime
            #print "     WHAM took: %f sec" % calctime
            
            #print startTF, round(startTF,1)        ## DEBUGGING
            temp = "%.1f" % round(startTF,1)
            outF = "free" + temp.split(".")[0] + temp.split(".")[1]
                
            a = np.loadtxt(outF,usecols=(1,))
            free_est[n,:len(a)] = a
            shutil.move(outF,"free_%d" % n)
            calc_iter = time.clock() - start_iter
            print "  iteration %d took: %f secs" % (n,calc_iter)

    free_std = np.std(free_est,axis=0)
    free_mean = np.mean(free_est,axis=0)
    #np.savetxt("free_bins.dat",np.loadtxt(outF,usecols=(0,)))
    np.savetxt("free_bins.dat",np.loadtxt("free_0",usecols=(0,)))
    np.savetxt("free_mean.dat",free_mean)
    np.savetxt("free_std.dat",free_std)     ## Use 2*free_std as error bars on free energy profile.
    os.chdir(cwd)

