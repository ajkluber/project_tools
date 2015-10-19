""" Compute Jacobian for matching a distance distribution 


Description:

    This module computes the jacobian of a distance distribution
such as measured with FRET.

note: as of now, only compute distances for FRET is updated

last updated: Justin Chen, May 05, 2015


"""

import numpy as np
import os
import time
import argparse

try:
    import mdtraj as md
except:
    pass

import model_builder as mdb

import scipy.stats as stats

from project_tools.parameter_fitting.util.util import *

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621
def_FRET_pairs = [[114,192]]
defspacing = 0.1 ## in nm


def find_sim_bins(savelocation, FRETr, fit_temp, residues=def_FRET_pairs, spacing=defspacing, weights=None):
    """find_sim_bins calculates and writes the simulation files """
    ##assumes nothing about where you are, but assumes the savelocation is specified correctly
    print "Calculating Simulation FRET bins"

    
    cwd = os.getcwd()
    
    if not os.path.isdir(savelocation):
        os.mkdir(savelocation)
     
    os.chdir(savelocation)
    
    #calcualte the parameters for binning
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = int(np.amin(FRETr)/spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*spacing,maxvalue*spacing)
    
    #if not weighted, set weights to ones
    if weights == None:
        weights = np.ones(np.shape(FRETr)[0])
    
    #actually histogram it
    print "***************************"
    print np.shape(FRETr)
    print np.shape(residues)
    print "***************************"
    hist, edges, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
    hist = hist/(np.sum(hist)*spacing)
    bincenters = 0.5 * (edges[1:] + edges[:-1])
    
    print "Making list of values:"
    #actually save it
    np.savetxt("simf_valuesT%d-P%d-%d.dat"%(fit_temp, residues[0], residues[1]),hist)
    np.savetxt("simf_edgesT%d-P%d-%d.dat"%(fit_temp, residues[0], residues[1]),edges)
    np.savetxt("simf-paramsT%d-P%d-%d.dat"%(fit_temp, residues[0], residues[1]),np.array([num_bins,minvalue*spacing,maxvalue*spacing,spacing]))
    
    os.chdir(cwd)
    print "Calculated bins for simulation data at a spacing of %.4f" % spacing
    
    return hist, slices

    
def fret_hist_calc(model, fitopts):
    fit_temp = fitopts["t_fit"]
    ##read trace file from 
    cwd = os.getcwd()
    subdir = model.name
    #load iteration number
    iteration = fitopts["iteration"]
    
    #load fret pairs and format correctly
    fret_pairs = fitopts["fret_pairs"]
    FRET_pairs = np.array(fret_pairs) - 1
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    FRETfile = "%s/FRET_hist.dat" % subdirec
    
    if not "fretdata" in fitopts:
        FRETtracefile = "%s/FRET_trace.dat" % cwd
    else:
        FRETtracefile = fitopts["fretdata"]
    
    
    for i in range(np.shape(FRET_pairs)[0]):
        residues = FRET_pairs[i,:]
        ftrace = np.loadtxt(FRETtracefile[0])
        edge_file = "%s/simf_edgesT%d-P%d-%d.dat"%(subdirec, fit_temp, residues[0], residues[1])
        edges = np.loadtxt(edge_file)
        hist, edges = np.histogram(ftrace, bins=edges, normed=True)
        bincenters = (edges[:-1] + edges[1:]) / 2
        datas = np.array([bincenters,hist])
        datas = np.transpose(datas)
        if i == 0:
            fret_total = datas
        else:
            fret_total = np.append(fret_total, datas, axis=0)
    np.savetxt(FRETfile, datas)
    
    print "Binned FRET_hist_calc Data"

def check_exp_data(FRETdata, bin_centers):
     #if correct within this marigin, then thinks its okay
     #Will check that the FRETdata centers and bin)centers are within 10^-6 of the spacing
    terms = np.shape(FRETdata)[0]
    i = 0
    
    spacing_FRET = FRETdata[1] - FRETdata[0]
    spacing_bin_centers = bin_centers[1] - bin_centers[0]
    
    min_difference = (min([spacing_FRET, spacing_bin_centers]))/1000000 #if correct within this marigin, then thinks its okay
    
    recalc = not np.shape(FRETdata)[0] == np.shape(bin_centers)[0]
    
    ##Verify that the bins line up 
    while (not recalc) and i<terms:
        if not (FRETdata[i] - bin_centers[i]) < min_difference:
            recalc = True
        i += 1
    return recalc

def add_error_log(note, fit_temp):
    errfile = "error_log-JC.txt"
    if not os.path.isfile(errfile):
        f = open("error_log-JC.txt","w")
        f.write("Error Log for This run\n\n")
        f.write("Global variables are:\n")
        f.write("Gas constant in kJ per mol = %d\n" % GAS_CONSTANT_KJ_MOL)
        f.write("pairs used are = " + str(def_FRET_pairs) + "\n")
        f.write("Temperature for Fitting used is T = %d\n" % fit_temp)
        f.write("Spacing in FRET pair distance used is = %d\n" %defspacing)
        f.write("\n")
        f.write(note)
        f.write("\n")
        f.close()
    else:
        f = f = open("error_log-JC.txt","a")
        f.write("\n")
        f.write(note)
        f.write("\n")
        f.close()
    print "ERROR: CHECK LOG \n %s" % note

        
def get_target_feature(model,fitopts):
    """ Get target features """
    fit_temp = fitopts["t_fit"]

        
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    FRETfile = "%s/FRET_hist.dat" % subdirec
        
    fret_hist_calc(model, fitopts)  
    FRETdata = np.loadtxt(FRETfile)  
    
    print "initial FRET data and bin_centers"
    print FRETdata
    
    target = FRETdata[:,1]
    target_err = target**0.5 ##for lack of a better way, take sqrt of bins for error estimate
    
    return target, target_err

def calculate_average_Jacobian(model,fitopts, FRET_pairs=def_FRET_pairs, spacing=defspacing ):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    if "t_fit" in fitopts:
        fit_temp = fitopts["t_fit"]
    else:
        raise IOError("Missing the fit_temperature, please specify in .ini file")
    
    if "fret_pairs" in fitopts:
        fret_pairs = fitopts["fret_pairs"]
        FRET_pairs = np.array(fret_pairs) - 1
        print "The FRET pairs are:"
        print FRET_pairs
    
    if "y_shift" in fitopts:
        y_shift = fitopts["y_shift"]   
    else:
        y_shift = 0.0
        fitopts["y_shift"] = 0.0        
 
    if "spacing" in fitopts:
        spacing = fitopts["spacing"]
    
    ##Define location of logical files    
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    traj_location = "%s/%d_0" % (sub, fit_temp)
    sim_location = "%s/fitting_%d" % (sub,iteration)
    ##define location of logical files
    
    os.chdir(traj_location)
    ## Get trajectory, state indicators, contact energy
    print "Working on calculating model's trajectory and contact info"
    traj,rij,qij = get_rij_Vp(model)

    ## Get simulation feature
    print "Now working on calculating the trajectories"
    beta = 1.0 / (GAS_CONSTANT_KJ_MOL*float(fit_temp))

    FRETr = md.compute_distances(traj,FRET_pairs, periodic=False)

    
    for i in range(np.shape(FRET_pairs)[0]):
        FRETr_use = FRETr[:,i] + y_shift
        print "Shifted simulated FRET-distance data by a y_shift = %f" % y_shift
        print FRETr_use
        sim_feature, sim_slices = find_sim_bins(sim_location, FRETr_use, fit_temp, residues=FRET_pairs[i,:], spacing=spacing, weights=None)
        print "Computing Jacobian and Simparams for the temperature %d, with spacing %f" % (fit_temp, spacing)
        Jacobian = compute_Jacobian_basic(qij,sim_feature*spacing, sim_slices, beta)
        Jacobian /= spacing        
        #store the sim_feature into a total array:
        if i == 0:
            sim_feature_all = sim_feature
            Jacobian_all = Jacobian
        else:
            sim_feature_all = np.append(sim_feature_all, sim_feature)
            Jacobian_all = np.append(Jacobian_all, Jacobian, axis=0)
    
    
    #save the temperature this was done in
    if not os.path.isdir("%s/newton"%sub):
        os.mkdir("%s/newton"%sub)
    f = open("%s/newton/temp-used-here.txt"%sub, "w")
    f.write(str(fit_temp))
    f.close()

    os.chdir(cwd)
    
    sim_feature_err = sim_feature_all ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian_all))
    
    return sim_feature_all, sim_feature_err, Jacobian_all, Jacobian_err

def compute_Jacobian_basic(qij, fr, sim_slices, beta, weights=None):
    """ Method for computing a Jacobian given only the rudimenary pieces necessary """
    ## qij is a NXM array containing the Qij values from the simulation
    ## fr is a RX1 array containing the normalized distributin f(r)
    ## Sim_slices is an NX1 array containing the bin_index+1 for the r matrix
    ## beta is the kbT for this particular Jacobian
    ## N = Number of frames, M = number of contacts to be fitted, R=number of bins of R data
    ## Note: assumes fr is already weighted!
    
    nbins = np.shape(fr)[0]
    (N_total_traj, npairs) = np.shape(qij)
    
    if weights == None:
        N_total_weight = N_total_traj
        weights = np.ones(N_total_traj)
    else:
        if not np.shape(weights)[0] == N_total_traj:
            raise IOError("Not every frame is weighted, aborting! Check to make sure weights is same length as trajectory")
        N_total_weight = np.sum(weights)
    
    Jacobian = np.zeros((nbins, npairs),float)
    for idx, bin_location in enumerate(sim_slices):
        Jacobian[bin_location-1, :] += qij[idx,:]*weights[idx]
    
    
    Jacobian /= N_total_weight
    Qavg = np.sum(Jacobian, axis=0)
    
    avg_matrix = np.dot(np.array([fr]).transpose(), np.array([Qavg]))
    print "The shape of these matrices are:"
    print np.shape(avg_matrix)
    print np.shape(Jacobian)
    Jacobian -= avg_matrix
    
    Jacobian *= (-1.0) * beta
    
    return Jacobian
    
    
if __name__ == "__main__":    
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='name.')
    parser.add_argument('--iteration', type=int, required=True, help='iteration.')
    args = parser.parse_args()
    
    name = args.name
    iteration= args.iteration

    pairs = np.loadtxt("%s/pairs.dat" % name,dtype=int)
    defaults = True
    model = mdb.models.SmogCalpha.SmogCalpha(name=name,pairs=pairs,defaults=defaults,iteration=iteration)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model)
