""" Compute Jacobian for matching a distance distribution 


Description:

    This module computes the jacobian of a distance distribution
such as measured with FRET.

note: as of now, only compute distances for FRET is updated

last updated: June 17, 2015


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

def calc_sim_bins(model, fitopts, framestep, residues=def_FRET_pairs, spacing=defspacing, weights=None):
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
    
    find_sim_bins(subdirec, FRETr[:,0], fit_temp, framestep, residues=residues, spacing=spacing, weights=weights)
    os.chdir(cwd)

def find_sim_bins(savelocation, FRETr, fit_temp, framestep, residues=def_FRET_pairs, spacing=defspacing, weights=None):
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
    hist, edges, slices = stats.binned_statistic(FRETr, weights, statistic="sum", range=[ran_size], bins=num_bins)
    hist = hist/(np.sum(hist)*spacing)
    bincenters = 0.5 * (edges[1:] + edges[:-1])
    
    # Get indices of transition bins
    F_indices, t_indices = get_transition_bins(slices, num_bins, framestep)
    
    print "Saving the bincenters:"
    #actually save it
    np.savetxt("simf_valuesT%d.dat"%fit_temp,hist)
    np.savetxt("simf_centers%d.dat"%fit_temp,bincenters)
    np.savetxt("simf-params%d.dat"%fit_temp,np.array([num_bins,minvalue*spacing,maxvalue*spacing,spacing]))
    
    os.chdir(cwd)
    print "Calculated bins for simulation data at a spacing of %.4f" % spacing
    
    return hist, F_indices, t_indices, num_bins
    
def get_transition_bins(slices, num_bins, framestep):
    # Get indices of transition bins
    # Returns linear indices of "transition bins," should match with numpy.ndarray.flatten() results
    # 'slices' has dimension Nx1, where N=number of frames
    
    # Get i and j slices, i=macrostate before transition, j=macrostate after transition
    # Subtract one to convert to zero-based indices
    F_indices = slices - 1
    state_i = F_indices[:-framestep]
    state_j = F_indices[framestep:]
    
    # Compute bin indices for transitions
    t_indices = state_i*num_bins + state_j
    
    return F_indices, t_indices

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
    
def tmatrix_exp_calc(model, fitopts, bin_size, ran_size, spacing):
    ##read trace file from 
    cwd = os.getcwd()
    subdir = model.name
    iteration = fitopts["iteration"]
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    TMfile = "%s/T_matrix_flat.dat" % subdirec
    
    ## MAY NEED TO MODIFY OPTIONS FOR T-MATRIX
    if not "TMdata" in fitopts:
        Tmatrixfile = "%s/T_matrix_exp.dat" % cwd
    else:
        Tmatrixfile = fitopts["TMdata"]
    
    Tmatrix = np.loadtxt(Tmatrixfile)

    #Extract entries matching range specified by ran_size
    lower_bin = int(np.around(ran_size[0]/spacing, 0))
    upper_bin = int(np.around(ran_size[1]/spacing, 0))
    T_matrix_small = Tmatrix[lower_bin:upper_bin, lower_bin:upper_bin]

    # Flatten and save
    T_matrix_flat = np.ndarray.flatten(T_matrix_small)
    np.savetxt(TMfile, T_matrix_flat)
    
    print "Extracted transition matrix"

def check_exp_data(TMdata, bin_centers):
    # Check to make sure dimensions match 
    ## ASSUMES EXPERIMENTAL T-MATRIX WAS PROPERLY BINNED, may need to change later if Tmatrix is rebinned every time
        
    spacing_FRET = TMdata[1] - TMdata[0]
    spacing_bin_centers = bin_centers[1] - bin_centers[0]
        
    recalc = not np.shape(TMdata)[0] == (np.shape(bin_centers)[0])**2
    
    return recalc
    
def find_FRET_bins(FRETr, framestep, spacing=defspacing):
    # Histograms the trace data into macrostate bins, analogous to find_sim_bins from compute_Jacobian script
    print np.shape(FRETr)
    weights = np.ones(np.shape(FRETr)[0])
    
    # Taken from find_sim_bins, included for consistency
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = int(np.amin(FRETr)/spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*spacing,maxvalue*spacing)
    
    hist, bins, slices = stats.binned_statistic(FRETr[:,0], weights, statistic="sum", range=[ran_size], bins=num_bins)
    
    # Call get_transition_bins to make transition bin indices
    F_indices, t_indices = get_transition_bins(slices, num_bins, framestep)
    
    return hist, F_indices, t_indices, num_bins
    
def get_T_Matrix(FRETr, framestep):    
    # Calculate flattened transition matrix for a given FRET trace
    # Based on fret_analysis/compute_transitions.py
    
    # Get FRET bins
    hist, F_indices, t_indices, num_bins = find_FRET_bins(FRETr, framestep)
    
    T_matrix = np.zeros((num_bins, num_bins))
    
    # Add ones to transition bins in square transition matrix
    for i in range(np.shape(F_indices)[0] - framestep):
        T_matrix[F_indices[i], F_indices[i+framestep]] += 1
    
    # Mask zeros to avoid divide-by-zero in normalization
    T_masked = np.ma.masked_where(T_matrix == 0, T_matrix)
    
    # Normalize each row
    for i in range(np.shape(T_matrix)[0]):
        T_masked[i,:] /= np.sum(T_masked[i,:])
        
    # Reshape to column vector
    T_matrix_flat = np.ndarray.flatten(T_masked)
    T_matrix_flat = np.transpose(T_matrix_flat)
    
    print np.shape(T_matrix_flat)
    
    return T_matrix_flat

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
    TMfile = "%s/T_matrix_flat.dat" % subdirec
    
    bin_centers = get_sim_centers(model, fitopts)
    bin_size, ran_size, spacing = get_sim_params(model, fitopts)
    ##Re-round to utilize the same range for ran_size as used in the simparams, in case floating point is round-off errored
    ran_size = (int(round(ran_size[0]/spacing))*spacing, int(round(ran_size[1]/spacing))*spacing)
    if not os.path.isfile(TMfile):
        tmatrix_exp_calc(model, fitopts, bin_size, ran_size, spacing)  
    TMdata = np.loadtxt(TMfile)  
    
    print "Initial Transition Matrix"

    print TMdata
    
    if check_exp_data(TMdata, bin_centers):
        print "Mismatched experimental data and simulated data. Attempting Re-binning"
        tmatrix_exp_calc(model, fitopts, bin_size, ran_size, spacing)
        TMdata = np.loadtxt(TMfile)
        if check_exp_data(TMdata, bin_centers):
            add_error_log("Catastrophic miscalculation of FRET bins", fit_temp)
            print "Found the TMdata to be:"
            print TMdata
    
    target = TMdata
        
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
    
    if "lag_step" in fitopts:
        framestep = fitopts["lag_step"]
    else:
        framestep = 1
        
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
    
    
    fr, F_indices, t_indices, num_bins = find_sim_bins(sim_location, FRETr[:,0], fit_temp, framestep, residues=FRET_pairs, spacing=spacing, weights=None)
    
    # FLATTENED T-MATRIX
    sim_feature = get_T_Matrix(FRETr, framestep)
    
    #calculate beta    
    beta = 1.0 / (GAS_CONSTANT_KJ_MOL*float(fit_temp))
    
    #save the temperature this was done in
    if not os.path.isdir("%s/newton"%sub):
        os.mkdir("%s/newton"%sub)
    f = open("%s/newton/temp-used-here.txt"%sub, "w")
    f.write(str(fit_temp))
    f.close()

    os.chdir(cwd)
    
    print "Computing Jacobian and T-matrix for the temperature %d, with spacing %f" % (fit_temp, spacing)
#    Tmatrix, Tmatrix_error = get_T_matrix(FRETr, ran_size, nbins) 
    Jacobian = compute_Jacobian_basic(qij, fr, framestep, F_indices, t_indices, beta)
    
    
    sim_feature_err = sim_feature ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian))
    
    return sim_feature, sim_feature_err, Jacobian, Jacobian_err

def compute_Jacobian_basic(qij, fr, framestep, F_indices, t_indices, beta, weights=None):
    """ Method for computing a Jacobian given only the rudimenary pieces necessary """
    ## qij is a NXM array containing the Qij values from the simulation
    ## fr is a RX1 array containing the normalized distributin f(r)
    ## Sim_slices is an NX1 array containing the bin_index+1 for the r matrix
    ## beta is the kbT for this particular Jacobian
    ## N = Number of frames, M = number of contacts to be fitted, R=number of bins of R data
    ## Note: assumes fr is already weighted!
    print "Starting compute_Jacobian_basic"
    # number of FRET bins
    nstates = np.shape(fr)[0]
    print qij[:,53]
    # number of TRANSITION bins
    nbins = nstates**2
        
    (N_total_traj, npairs) = np.shape(qij)
    
    qi = np.zeros((nstates, npairs))
    qi_count = np.zeros(nstates)
        
    if weights == None:
        N_total_weight = N_total_traj
        weights = np.ones(N_total_traj)
    else:
        if not np.shape(weights)[0] == N_total_traj:
            raise IOError("Not every frame is weighted, aborting! Check to make sure weights is same length as trajectory")
        N_total_weight = np.sum(weights)

    # Q values for all frames starting in a particular bin
    for idx, F_bin_location in enumerate(F_indices[:-framestep]):
        qi[F_bin_location,:] += (qij[idx,:] + qij[idx+framestep,:])
        qi_count[F_bin_location] += 1
    
    #normalize the average value for the pair sum starting in state i
    for i in range(np.shape(qi)[0]):
        if qi_count[i] != [0]:
            qi[i,:] /= float(qi_count[i])

    Jacobian = np.zeros((nbins, npairs))
    for idx, t_bin_location in enumerate(t_indices):
        # Add q values for specific transition
        Jacobian[t_bin_location, :] += (qij[idx,:] + qij[idx+framestep,:])    
    print "*********"
    print Jacobian[10,53]
    # Normalize by i bin count
    for i in range(np.shape(Jacobian)[0]):
        state_i_idx = int(np.floor(i/nstates))
        if qi_count[state_i_idx] != 0:
            Jacobian[i,:] /= (qi_count[state_i_idx])
    print "*********"
    print Jacobian[10,53]
    for idx, t_bin_location in enumerate(t_indices):
        # Index for q value of all transitions starting at state i
        state_i_idx = int(np.floor(t_bin_location/nstates))
        # Subtract q values for all starting at state i
        Jacobian[t_bin_location, :] -= qi[state_i_idx,:]
    
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
