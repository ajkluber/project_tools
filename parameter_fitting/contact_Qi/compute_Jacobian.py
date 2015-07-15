""" Compute Jacobian for contact probability function


Description:

    This module computes the jacobian of the contact probability function.

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

#import project_tools.parameter_fitting.util.util as util
from project_tools.parameter_fitting.util.util import *

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def get_target_feature(model):
    """ Get target features """
    name = model.subdir
    iteration = model.iteration
    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    
    ## To Do:
    ## - Check if a target set of contact probabilities is given
    ##   else construct target as <Q_i^TS> = Q^TS (uniform TS).

    ## ---- For future
    ## Format for target_Qi.dat
    ## - three columns:
    ##  <res_a>  <res_b>  <Q_ab^TS>
    ## computes contact as within native contact distance.
    ## ---- 
    
    ## ---- For now 
    ## Just a column with the desired contact probability in the TS

    if os.path.exists("%s/target_Qi.dat" % name):
        target =  np.loadtxt("%s/target_Qi.dat" % name)
        target_err =  np.loadtxt("%s/target_Qi_err.dat" % name)
    else:
        ## Compute the average Q of the TS: Average of the endpoints.
        os.chdir("%s" % sub)
        #bounds, state_labels = util.get_state_bounds()
        bounds, state_labels = get_state_bounds()
        Q_TS = 0.5*(bounds[2] + bounds[3])/float(model.n_contacts)
        target = Q_TS*np.ones(model.n_contacts,float)
        target_err = 0.05*np.ones(model.n_contacts,float)
        os.chdir(cwd)

    return target, target_err

def calculate_average_Jacobian(model):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    
    name = model.subdir
    iteration = model.iteration

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("long_temps_last","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 

    #bounds, state_labels = util.get_state_bounds()
    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in iteration subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. 
    sim_feature_all = []
    Jacobian_all = []
    lasttime = time.time()
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for iteration_%d/%s" % (model.Mut_iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,beta,bounds)
        sim_feature_all.append(sim_feature)
        Jacobian_all.append(Jacobian)
        os.chdir("..")

        thistime = time.time()
        timediff = thistime - lasttime
        lasttime = thistime
        print "  calculation took %.2f seconds = %.2f minutes" % (timediff,timediff/60.)

    sim_feature_all = np.array(sim_feature_all)
    Jacobian_all = np.array(Jacobian_all)

    ## Take avg. and use standard deviation as error bars.
    sim_feature_avg = sum(sim_feature_all)/float(len(directories))
    sim_feature_err = np.std(sim_feature_all,axis=0)
    Jacobian_avg = sum(Jacobian_all)/float(len(directories))
    Jacobian_err = np.std(Jacobian_all,axis=0)

    os.chdir(cwd)

    return sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err

def compute_Jacobian_for_directory(model,beta,bounds):
    """ Calculates the feature vector (ddG's) and Jacobian for one directory """

    ## Get trajectory, state indicators, contact energy
    Q = np.loadtxt("Q.dat")
    #traj,rij,Vp = util.get_rij_Vp(model)
    #U,TS,N,Uframes,TSframes,Nframes = util.get_state_indicators(Q,bounds)
    traj,rij,Vp = get_rij_Vp(model)
    U,TS,N,Uframes,TSframes,Nframes = get_state_indicators(Q,bounds)

    ## Get Qi
    Qi = np.loadtxt("qimap.dat",dtype=float)
    sim_feature = sum(Qi[TS,:])/TSframes

    ## Initialize Jacobian
    #Jacobian = np.zeros((model.n_contacts,model.n_contacts),float)
    Jacobian = np.zeros((model.n_contacts,model.n_model_param),float)

    avg_Vp = sum(Vp[TS,:])/TSframes
    ## Compute rows of the Jacobian which are correlation functions of 
    ## contact formation with contact energy.
    for i in range(model.n_contacts):
        if (i % 10) == 0:
            print "    row %d out of %d" % (i+1,model.n_contacts)
        avg_Qi = sum(Qi[TS,i])/TSframes 
        for j in range(model.n_model_param):
            avg_QiVp = sum(Qi[TS,i]*Vp[TS,j])/TSframes
            Jacobian[i,j] = -beta*(avg_QiVp - avg_Qi*avg_Vp[j])

    return sim_feature, Jacobian

def calculate_average_Jacobian(model,fitopts):
    """ Calculate ddG's and Jacobian over replicas"""
    
    name = model.name
    iteration = fitopts['iteration']

    # Get mutations and fraction of native pairs deleted for each mutation.
    os.chdir("%s/iteration_%d" % (name,iteration))


    Tlist = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 
    Tlist = Tlist[0]
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(Tlist[0].split("_")[0]))
    
    # Get pairwise distances from trajectories
    trajfiles = [ "%s/traj.xtc" % x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 
    native = "%s/Native.pdb" % Tlist[0]
    rij = get_rij(model,trajfiles,native)
    contact_distance = np.array([ model.pairwise_other_parameters[x][0] for x in range(1,model.n_pairs,2) ])

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_pairs]

    # Found speedups by calculating quantities on per state basis
    U,TS,N,Uframes,TSframes,Nframes = concatenate_state_indicators(Tlist,bounds,coord="Q.dat")

    # Average dimensionless potential energy for each state
    Vp_U   = get_Vp_for_state(model,rij,U,Uframes)
    Vp_TS  = get_Vp_for_state(model,rij,TS,TSframes)
    Vp_N   = get_Vp_for_state(model,rij,N,Nframes)
    sumVp_U = np.mean(Vp_U,axis=0); sumVp_TS = np.mean(Vp_TS,axis=0); sumVp_N = np.mean(Vp_N,axis=0) 

if __name__ == "__main__":    
    """Calculate the  """
    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='Directory.')
    #parser.add_argument('--iteration', type=int, required=True, help='Iteration.')
    args = parser.parse_args()

    name = args.name
    
    model, fitopts = mdb.inputs.load_model(name)
    iteration = fitopts['iteration']
    model.fitting_params = np.arange(1,model.n_pairs,2)
    model.n_fitting_params = len(model.fitting_params)
    
    wantpairs = model.pairs[1::2] - 1

    # Get mutations and fraction of native pairs deleted for each mutation.
    os.chdir("%s/iteration_%d" % (name,iteration))

    Tlist = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 
    Tlist = [ Tlist[0] ]
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(Tlist[0].split("_")[0]))
    
    # Get pairwise distances from trajectories
    trajfiles = [ "%s/traj.xtc" % x.rstrip("\n") for x in Tlist ] 
    native = "%s/Native.pdb" % Tlist[0]
    n_residues = len(open(native,"r").readlines()) - 1
    rij = get_rij(model,trajfiles,native)
    contact_r0 = np.array([ model.pairwise_other_parameters[x][0] for x in model.fitting_params ])
    
    # Get state bounds
    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_pairs]

    # Found speedups by calculating quantities on per state basis
    U,TS,N,Uframes,TSframes,Nframes = concatenate_state_indicators(Tlist,bounds,coord="Q.dat")

    # Average potential energy for each state
    Vp_U   = beta*get_Vp_for_state(model,rij,U,Uframes)
    Vp_TS  = beta*get_Vp_for_state(model,rij,TS,TSframes)
    Vp_N   = beta*get_Vp_for_state(model,rij,N,Nframes)
    sumVp_U = np.mean(Vp_U,axis=0); sumVp_TS = np.mean(Vp_TS,axis=0); sumVp_N = np.mean(Vp_N,axis=0) 

    # Average of contact function for each state
    rij_nat = rij[:,model.fitting_params] 
    Q_U  = (rij_nat[U,:]  < (contact_r0 + 0.1)).astype(int)
    Q_TS = (rij_nat[TS,:] < (contact_r0 + 0.1)).astype(int)
    Q_N  = (rij_nat[N,:]  < (contact_r0 + 0.1)).astype(int)
    sumQ_U = np.mean(Q_U,axis=0); sumQ_TS = np.mean(Q_TS,axis=0); sumQ_N = np.mean(Q_N,axis=0)

    # Jacobian
    J_U  = np.dot(Q_U.T,Vp_U)   - np.outer(sumQ_U,sumVp_U)
    J_TS = np.dot(Q_TS.T,Vp_TS) - np.outer(sumQ_TS,sumVp_TS)
    J_N  = np.dot(Q_N.T,Vp_N)   - np.outer(sumQ_N,sumVp_N)

    # Hessian
    H_U  = np.dot(J_U.T,J_U)
    H_TS = np.dot(J_TS.T,J_TS)
    H_N  = np.dot(J_N.T,J_N)

    vals_U, vecs_U = np.linalg.eig(H_U)    
    vals_TS, vecs_TS = np.linalg.eig(H_TS)    
    vals_N, vecs_N = np.linalg.eig(H_N)    

    C_U = np.zeros((n_residues,n_residues))
    for i in range(model.n_native_pairs):
        C_U[model.native_pairs[i][1] - 1, model.native_pairs[i][0] - 1] = vecs_U[i,0]
    C_TS = np.zeros((n_residues,n_residues))
    for i in range(model.n_native_pairs):
        C_TS[model.native_pairs[i][1] - 1, model.native_pairs[i][0] - 1] = vecs_TS[i,0]
    C_N = np.zeros((n_residues,n_residues))
    for i in range(model.n_native_pairs):
        C_N[model.native_pairs[i][1] - 1, model.native_pairs[i][0] - 1] = vecs_N[i,0]

    os.chdir("../..")
