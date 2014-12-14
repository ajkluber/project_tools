""" Compute Jacobian for matching a distance distribution 


Description:

    This module computes the jacobian of a distance distribution
such as measured with FRET.

note: as of now, only compute distances for FRET is updated

last updated: Justin Chen, December 03, 2014


"""

import numpy as np
import os
import time
import argparse

import mdtraj as md

import model_builder as mdb

from project_tools.parameter_fitting.util.util import *

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621
FRET_pairs = [[114,192]]

def return_temp_file_value():
    if os.path.isfile("fitting_temperature.txt"):
        print "Setting new temperature based on fitting_temperature.txt"
        return np.loadtxt("fitting_temperature.txt")
    else:
        print "No fitting_temperature specified, using default T=0"
        return 0
        
def_temp = return_temp_file_value()
defspacing = 0.1 ## in nm


def calc_sim_bins(model,residues=FRET_pairs,fit_temp=def_temp,spacing=defspacing):
    """calc_sim_bins calculates and writes the simulation files """
    ##assumes you are in the folder containing the model subdir
    print "Calculating Simulation FRET bins"
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subtemp = "%s/%d_%d" % (sub,fit_temp,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    
    if not os.path.isdir(subdirec):
        os.mkdir(subdirec)
     
    os.chdir(subtemp)
    traj = md.load("traj.xtc",top="Native.pdb")
    
    os.chdir(subdirec)
    FRETr = md.compute_distances(traj,residues, periodic=False)
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = int(np.amin(FRETr)/spacing)
    num_bins = maxvalue - minvalue
    ran_size = (minvalue*spacing,maxvalue*spacing)
    hist, edges = np.histogram(FRETr,num_bins,ran_size)
    hist = hist/(np.sum(hist)*spacing)
    np.savetxt("simf_valuesT%d.dat"%fit_temp,hist)
    edges = edges + (0.5*spacing)
    bincenters = edges[:-1]
    bincenters = np.arange((minvalue*spacing)+(spacing/2.0), maxvalue*spacing, spacing)
    np.savetxt("simf_centers%d.dat"%fit_temp,bincenters)
    np.savetxt("simf-params%d.dat"%fit_temp,np.array([num_bins,minvalue*spacing,maxvalue*spacing,spacing]))
    
    os.chdir(cwd)
    print "Calculated bins for simulation data at a spacing of %.4f" % spacing

def get_sim_params(model,fit_temp=def_temp):
    ##assumes you are in the folder containing the model subdir
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    parmfile = "%s/simf-params%d.dat" % (subdirec,fit_temp)
    
    if not os.path.isfile(parmfile):
        calc_sim_bins(model,fit_temp=fit_temp)  
    parms = np.loadtxt(parmfile)
    num_bins = parms[0]
    ran_size = (parms[1],parms[2])
    spacing = parms[3]
    return num_bins, ran_size, spacing

def get_sim_centers(model,fit_temp=def_temp): 
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    if not os.path.isfile(simfile):
        calc_sim_bins(model,fit_temp=fit_temp)
    return np.loadtxt(simfile)

def get_sim_array(model,fit_temp=def_temp):
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    simfilef = "%s/simf_valuesT%d.dat" % (subdirec,fit_temp)
    print "Fit temp is : ", fit_temp
    if not os.path.isfile(simfile):
        calc_sim_bins(model,fit_temp=fit_temp)
    centers = np.loadtxt(simfile)
    values = np.loadtxt(simfilef)
    simarray = np.array([centers,values])
    simarray = np.transpose(simarray)
    
    os.chdir(cwd)
    return simarray
    
def fret_hist_calc(model, bin_size, ran_size, spacing):
    ##read trace file from 
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    FRETfile = "%s/FRET_hist.dat" % subdirec
    FRETtracefile = "%s/FRET_trace.dat" % cwd
    
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
    terms = np.shape(FRETdata)[0]
    i = 0
    recalc = False
    ##Verify that the bins line up 
    while (not recalc) and i<terms:
        if not FRETdata[i] == bin_centers[i]:
            recalc = True
        i += 1
    return recalc

def add_error_log(note):
    errfile = "error_log-JC.txt"
    if not os.path.isfile(errfile):
        f = open("error_log-JC.txt","w")
        f.write("Error Log for This run\n\n")
        f.write("Global variables are:\n")
        f.write("Gas constant in kJ per mol = %d\n" % GAS_CONSTANT_KJ_MOL)
        f.write("pairs used are = " + str(FRET_pairs) + "\n")
        f.write("Temperature for Fitting used is T = %d\n" % def_temp)
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
    print "ERROR: CHECK LOG \n"

        
def get_target_feature(model,fit_temp=def_temp):
    """ Get target features """
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    subdirec = "%s/fitting_%d" % (sub,iteration)
    simfile = "%s/simf_centers%d.dat" % (subdirec,fit_temp)
    FRETfile = "%s/FRET_hist.dat" % subdirec
    FRETtracefile = "%s/FRET_trace.dat" % sub
    
    bin_centers = get_sim_centers(model)
    bin_size, ran_size, spacing = get_sim_params(model)
   
    if not os.path.isfile(FRETfile):
        fret_hist_calc(model, bin_size, ran_size, spacing)  
    FRETdata = np.loadtxt(FRETfile)  
    
    print FRETdata[:,0]
    print bin_centers
    if check_exp_data(FRETdata[:,0], bin_centers):
        print "Mismatched experimental data and simulated data. Attempting Re-binning"
        fret_hist_calc(model, bin_size, ran_size, spacing)
        FRETdata = np.loadtxt(FRETfile)
        if check_exp_data(FRETdata[:,0], bin_centers):
            add_error_log("Catastrophic miscalculation of FRET bins")
    
    target = FRETdata[:,1]
    target_err = target**0.5 ##for lack of a better way, take sqrt of bins for error estimate
    
    return target, target_err

def calculate_average_Jacobian(model,residues=FRET_pairs,fit_temp=def_temp):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    cwd = os.getcwd()
    subdir = model.subdir
    iteration = model.iteration
    sub = "%s/%s/iteration_%d" % (cwd,subdir,iteration)
    os.chdir(sub)
    os.chdir("%d_%d" % (fit_temp,iteration))
    beta = 1.0 * (GAS_CONSTANT_KJ_MOL*float(fit_temp))
    
    
    #save the temperature this was done in
    if not os.path.isdir("%s/newton"%sub):
        os.mkdir("%s/newton"%sub)
    f = open("%s/newton/temp-used-here.txt"%sub, "w")
    f.write(str(fit_temp))
    f.close()
    
    print "Computing Jacobian and Simparams for the temperature %d, with spacing %f" % (fit_temp, defspacing)
    Jacobian, simparams = compute_Jacobian_for_directory(model,beta,FRET_pairs,defspacing)
    
    os.chdir(cwd)
    if not np.array_equal(simparams, get_sim_array(model,fit_temp)):
        add_error_log("Catastrophic error. Simparams from compute_Jacobian_for_directory is different from get_sim_params")
    
    sim_feature = simparams[:,1]
    sim_feature_err = sim_feature ** 0.5
    Jacobian_err = np.zeros(np.shape(Jacobian))
    
    return sim_feature, sim_feature_err, Jacobian, Jacobian_err

def compute_Jacobian_for_directory(model,beta,residues,spacing):
    """ Calculates the Jacobian for one directory """

    ## Get trajectory, state indicators, contact energy
    print "Working on calculating model's trajectory and contact info"
    traj,rij,qij = get_rij_Vp(model)

    ## Get simulation feature
    print "Now working on calculating the trajectories"
    FRETr = md.compute_distances(traj,residues, periodic=False)
   
    ## Initialize Jacobian
    print "Now preparing the Jacobians"
    
    ##Find the fret Trajectories largest and smallest value
    maxvalue = int(np.amax(FRETr)/spacing) + 1
    minvalue = int(np.amin(FRETr)/spacing)
    total_traj = np.shape(rij)[0]
    num_bins = maxvalue - minvalue
    
    Jacobian = np.zeros((num_bins,model.n_contacts),float)
    Fr = np.zeros(num_bins, float)  ##Number of counts of f for that r
    Qr = np.zeros((num_bins,model.n_contacts), float)
    Qcount = np.zeros(num_bins, int)  ##Total value of Qs for that r
    Qavg = np.zeros(model.n_contacts)   ##Average value of Q over whole trajectory 
    
    ## Jacobian 1st cord = F_r, probability-dist at different r
    ## Jacobian 2nd cord = Derivative wrt native contact epsilon
    ## Compute rows of the Jacobian which are correlation functions of 
    ## contact formation with contact energy.
    count = 0
    print "Warning, starting Jacobian Calculation. ONLY configured for one pair"
    for i in FRETr:
        FRETbin = int(i/spacing) - minvalue
        Qcount[FRETbin] += 1
        Qr[FRETbin,:] += qij[count,:]
        Qavg += qij[count,:]
        Fr[FRETbin] += 1
        count += 1
    
    ##error checking
    if np.sum(Qcount) > total_traj:
        print "ERROR: Number of binned frames greater than total number of frames"
    elif np.sum(Qcount) < total_traj:
        print "ERROR: Number of binned frames less than total number of frames"

    print "Beginning Normalization and Averaging, Standby"
    count = 0
    for i in Qcount:
        if not i == 0:
            Qr[count] /= i
            count += 1
    
    Qavg /= total_traj
    Fr /= np.sum(Fr)*spacing
    
    print "Calculating the Final Jacobian, Please Standby ... "
    
    
    Jacobian = np.zeros((num_bins,model.n_contacts), float)
    Jacobian += Qr
    
    countr = 0
    for i in Fr:
        counte = 0
        for j in Qavg:
            Jacobian[countr,counte] += i*j
            counte += 1
        countr =+ 1
    Jacobian *= -1 * beta
    
    xFRET = np.arange((minvalue*spacing)+(spacing/2.0), maxvalue*spacing, spacing)
    simparams = np.array([xFRET,Fr])
    simparams = np.transpose(simparams)
    print "Finished computing the Jacobian and Simparams"
    return Jacobian, simparams

if __name__ == "__main__":    
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='name.')
    parser.add_argument('--iteration', type=int, required=True, help='iteration.')
    args = parser.parse_args()
    
    name = args.name
    iteration= args.iteration

    contacts = np.loadtxt("%s/contacts.dat" % name,dtype=int)
    pdb = "%s.pdb" % name
    defaults = True
    model = mdb.models.SmogCalpha.SmogCalpha(pdb=pdb,contacts=contacts,defaults=defaults,iteration=iteration)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model)
