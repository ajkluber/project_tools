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

def calc_sim_bins(model, fitopts, residues=def_FRET_pairs, spacing=defspacing, weights=None):
    """calc_sim_bins does the actual calculation, with minimal assumptions of directory location """
    ##assumes you are in the folder containing the model subdir
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
        y_shift = 0
        fitopts["y_shift"] = 0        
 
    if "spacing" in fitopts:
        spacing = fitopts["spacing"]
        
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    fit_temp = fitopts["t_fit"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subtemp = "%s/%d_0" % (sub,fit_temp)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    
    if not os.path.isdir(subdirec):
        os.mkdir(subdirec)
     
    os.chdir(subtemp)
    traj = md.load("traj.xtc",top="Native.pdb")
    FRETr = md.compute_distances(traj,residues, periodic=False)
    FRETr += y_shift 
    
    print "FRETr after md.compute is:"
    print FRETr
    print np.shape(FRETr)
    
    find_sim_bins(subdirec, FRETr[:,0], fit_temp, residues=residues, spacing=spacing, weights=weights)
    os.chdir(cwd)

def find_sim_bins(savelocation, FRETr, fit_temp, residues=def_FRET_pairs, spacing=defspacing, weights=None):
    """find_sim_bins calculates and writes the simulation files """
    ##assumes nothing about where you are, but assumes the savelocation is specified correctly
    print "Calculating Simulation FRET bins"
    #savelocation should be in "cwd/subdir/iteration_number/fitting_number
    
    ##debugging
    print "save location is in:"
    print savelocation
    print "FRETr is:"
    print FRETr
    print "shape of FRETr is:"
    print np.shape(FRETr)
    print "fit temp is:"
    print fit_temp
    print "residues are:"
    print residues
    print "spacing is:"
    print spacing
    print "weights is"
    print weights
    ##debugging
    
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
    print "***************************"
    hist, edges, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
    hist = hist/(np.sum(hist)*spacing)
    bincenters = 0.5 * (edges[1:] + edges[:-1])
    
    print "Saving the bincenters:"
    #actually save it
    np.savetxt("simf_valuesT%d.dat"%fit_temp,hist)
    np.savetxt("simf_centers%d.dat"%fit_temp,bincenters)
    np.savetxt("simf-params%d.dat"%fit_temp,np.array([num_bins,minvalue*spacing,maxvalue*spacing,spacing]))
    
    os.chdir(cwd)
    print "Calculated bins for simulation data at a spacing of %.4f" % spacing
    
    return hist, slices

def get_sim_params(model,fitopts):
    ##assumes you are in the folder containing the model subdir
    fit_temp = fitopts["t_fit"]
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    parmfile = "%s/simf-params%d.dat" % (subdirec,fit_temp)
    
    if not os.path.isfile(parmfile):
        calc_sim_bins(model, fitopts)  
    parms = np.loadtxt(parmfile)
    num_bins, ran_size, spacing = analyze_sim_params(parms)
    return num_bins, ran_size, spacing

def analyze_sim_params(parms):
    num_bins = parms[0]
    ran_size = (parms[1],parms[2])
    spacing = parms[3]
    return num_bins, ran_size, spacing

def get_sim_centers(model,fitopts): 
    fit_temp = fitopts["t_fit"]
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    if not os.path.isfile(simfile):
        find_sim_bins(model,fitopts)
    return np.loadtxt(simfile)

def get_sim_array(model,fitopts):
    fit_temp = fitopts["t_fit"]
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    simfilef = "%s/simf_valuesT%d.dat" % (subdirec,fit_temp)
    print "Fit temp is : ", fit_temp
    if not os.path.isfile(simfile):
        find_sim_bins(model,fitopts)
    centers = np.loadtxt(simfile)
    values = np.loadtxt(simfilef)
    simarray = np.array([centers,values])
    simarray = np.transpose(simarray)
    
    os.chdir(cwd)
    return simarray
    
def fret_hist_calc(model, fitopts, bin_size, ran_size, spacing):
    ##read trace file from 
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    FRETfile = "%s/FRET_hist.dat" % subdirec
    if not "fretdata" in fitopts:
        FRETtracefile = "%s/FRET_trace.dat" % cwd
    else:
        FRETtracefile = fitopts["fretdata"][0]
    
    ftrace = np.loadtxt(FRETtracefile)
    hist, edges = np.histogram(ftrace,bin_size,ran_size)
    #normalize and find bin centers
    hist = hist/(np.shape(ftrace)[0]*spacing)
    edges = edges + (0.5*spacing)
    edges = edges[:-1]
    bincenters = np.arange(ran_size[0]+(spacing/2.0), ran_size[1], spacing)
    datas = np.array([bincenters,hist])
    datas = np.transpose(datas)
    np.savetxt(FRETfile, datas)
    
    print "Rebinned FRET_hist_calc Data"

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
    
    bin_centers = get_sim_centers(model, fitopts)
    bin_size, ran_size, spacing = get_sim_params(model, fitopts)
    ##Re-round to utilize the same range for ran_size as used in the simparams, in case floating point is round-off errored
    ran_size = (int(round(ran_size[0]/spacing))*spacing, int(round(ran_size[1]/spacing))*spacing)
    if not os.path.isfile(FRETfile):
        fret_hist_calc(model, fitopts, bin_size, ran_size, spacing)  
    FRETdata = np.loadtxt(FRETfile)  
    
    print "initial FRET data and bin_centers"

    print FRETdata[:,0]
    print bin_centers
    
    if check_exp_data(FRETdata[:,0], bin_centers):
        print "Mismatched experimental data and simulated data. Attempting Re-binning"
        fret_hist_calc(model, fitopts, bin_size, ran_size, spacing)
        FRETdata = np.loadtxt(FRETfile)
        if check_exp_data(FRETdata[:,0], bin_centers):
            add_error_log("Catastrophic miscalculation of FRET bins", fit_temp)
            print "Found the FRETdata and bin_centers to be:"
            print FRETdata[:,0]
            print bin_centers
    
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
    FRETr = md.compute_distances(traj,FRET_pairs, periodic=False)
    print FRETr
    FRETr += y_shift
    print "Shifted simulated FRET-distance data by a y_shift = %f" % y_shift
    print FRETr
    ##WARNING: CONFIGURED FOR ONLY ONE PAIR OF FRET PROBES
    ## To configure for more, modify such that you iterate over each pair in FRETr, calcualte the sim_feature, sim_slices
    ## The nrun the Jacobian on just that slice
    ## then rerun sim_feature, sim_slices using next pair
    ## append at the end
    sim_feature, sim_slices = find_sim_bins(sim_location, FRETr[:,0], fit_temp, residues=FRET_pairs, spacing=spacing, weights=None)
    
    beta = 1.0 / (GAS_CONSTANT_KJ_MOL*float(fit_temp))
    
    
    #save the temperature this was done in
    if not os.path.isdir("%s/newton"%sub):
        os.mkdir("%s/newton"%sub)
    f = open("%s/newton/temp-used-here.txt"%sub, "w")
    f.write(str(fit_temp))
    f.close()

    os.chdir(cwd)
    
    print "Computing Jacobian and Simparams for the temperature %d, with spacing %f" % (fit_temp, spacing)
    Jacobian = compute_Jacobian_basic(qij,sim_feature*spacing, sim_slices, beta)
    Jacobian /= spacing        
    sim_feature_err = sim_feature ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian))
    
    return sim_feature, sim_feature_err, Jacobian, Jacobian_err

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
