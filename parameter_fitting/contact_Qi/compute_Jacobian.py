""" Compute Jacobian for contact probability function


Description:

    This module computes the jacobian of the contact probability function.

DEPRECATED
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

if __name__ == "__main__":    
    
    name = "1FMK"
    #iteration = 1

    model = mdb.check_inputs.load_model(name)
    #iteration = model.iteration
    iteration = 0

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
    #for n in range(len(directories)):
    lasttime = time.time()
    for n in [0]:
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for iteration_%d/%s" % (model.iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,beta,bounds)

        thistime = time.time()
        timediff = thistime - lasttime
        lasttime = thistime
        print "  calculation took %.2f seconds = %.2f minutes" % (timediff,timediff/60.)

    np.savetxt("test_J.dat",Jacobian)
    os.chdir(cwd)


    """
    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='Directory.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration.')
    args = parser.parse_args()
    name = args.name
    iteration = args.iteration
    model = mdb.check_inputs.load_model(name) 
    model.iteration = iteration
    target_feature, target_feature_err = get_target_feature(model)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model)

    if not os.path.exists("%s/iteration_%d/contact_Qi" % (name,iteration)):
        os.mkdir("%s/iteration_%d/contact_Qi" % (name,iteration))

    np.savetxt("%s/iteration_%d/contact_Qi/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/iteration_%d/contact_Qi/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/iteration_%d/contact_Qi/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/iteration_%d/contact_Qi/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/iteration_%d/contact_Qi/sim_feature.dat" % (name,iteration), sim_feature_avg)
    np.savetxt("%s/iteration_%d/contact_Qi/sim_feature_err.dat" % (name,iteration), sim_feature_err)
    np.savetxt("%s/iteration_%d/contact_Qi/Jacobian.dat" % (name,iteration), Jacobian_avg)
    np.savetxt("%s/iteration_%d/contact_Qi/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)
    """
