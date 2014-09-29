""" Compute Jacobian for contact probability function


Description:

    This module computes the jacobian of the contact probability function.

"""

import numpy as np
import os
import time
import argparse

import mdtraj as md

import model_builder as mdb

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def get_state_bounds():
    """ Bounds for each state. Bounds are bin edges along Q. """
    if os.path.exists("state_bounds.txt"):
        statefile = open("state_bounds.txt","r").readlines()
    else:
        print "ERROR!"
        print "  Please create state_bounds.txt"
        print "  With the boundaries of each state along Q"
        print "  Exiting"
        raise SystemExit
    
    state_bounds = []
    state_labels = []
    for line in statefile:
        info = line.split()
        state_labels.append(info[0])
        state_bounds.append(float(info[1]))
        state_bounds.append(float(info[2]))
    
    return state_bounds,state_labels

def get_states_Vij(model,bounds):
    """ Load trajectory, state indicators, and contact energy """

    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape))
    Q = np.loadtxt("Q.dat") ## To Do: Generalize to different reaction coordinates

    state_indicator = np.zeros(len(Q),int)
    ## Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (Q > bounds[state_num]).astype(int)*(Q <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        num_not_assign = sum((state_indicator == 0).astype(int))
        print "  Warning! %d frames were not assigned out of %d total frames!" % (num_not_assign,len(Q))
    ## Boolean arrays that indicate which state each frame is in.
    ## States are defined by their boundaries along coordinate Q.
    U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
    TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
    N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))

    Vij = model.calculate_contact_potential(rij,"all")

    return traj,U,TS,N,Uframes,TSframes,Nframes,Vij

def get_target_feature(model):
    """ Get target features """
    name = model.subdir
    iteration = model.Mut_iteration
    cwd = os.getcwd()
    sub = "%s/%s/Mut_%d" % (cwd,name,iteration)
    
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
        bounds, state_labels = get_state_bounds()
        Q_TS = 0.5*(bounds[2] + bounds[3])/float(model.n_contacts)
        target = Q_TS*np.ones(model.n_contacts,float)
        target_err = 0.05*np.ones(model.n_contacts,float)
        os.chdir(cwd)

    return target, target_err

def calculate_average_Jacobian(model):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    
    name = model.subdir
    iteration = model.Mut_iteration

    cwd = os.getcwd()
    sub = "%s/%s/Mut_%d" % (cwd,name,iteration)
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in Mut subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. 
    sim_feature_all = []
    Jacobian_all = []
    lasttime = time.time()
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for Mut_%d/%s" % (model.Mut_iteration,dir)
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
    traj,U,TS,N,Uframes,TSframes,Nframes,Vij = get_states_Vij(model,bounds)

    ## Get Qi
    Qi = np.loadtxt("qimap.dat",dtype=float)
    sim_feature = sum(Qi[TS,:])/TSframes

    ## Initialize Jacobian
    Jacobian = np.zeros((model.n_contacts,model.n_contacts),float)

    ## Compute rows of the Jacobian which are correlation functions of 
    ## contact formation with contact energy.
    for i in range(model.n_contacts):
        if (i % 10) == 0:
            print "    row %d out of %d" % (i+1,model.n_contacts)

        avg_Qi = sum(Qi[TS,i])/TSframes 
        avg_Vij = sum(Vij[TS,:])/TSframes
        avg_QiVij = sum((Qi[TS,i]*Vij[TS,:].T).T)/TSframes

        Jacobian[i,:] = -beta*(avg_QiVij - avg_Qi*avg_Vij)

    return sim_feature, Jacobian

if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='Directory.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    model = mdb.check_inputs.load_model(name) 
    model.Mut_iteration = iteration
    
    target_feature, target_feature_err = get_target_feature(model)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model)

    if not os.path.exists("%s/Mut_%d/contact_Qi" % (name,iteration)):
        os.mkdir("%s/Mut_%d/contact_Qi" % (name,iteration))

    np.savetxt("%s/Mut_%d/contact_Qi/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/Mut_%d/contact_Qi/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/Mut_%d/contact_Qi/target_feature.dat" % (name,iteration), target_feature)
    np.savetxt("%s/Mut_%d/contact_Qi/target_feature_err.dat" % (name,iteration), target_feature_err)
    np.savetxt("%s/Mut_%d/contact_Qi/sim_feature.dat" % (name,iteration), sim_feature_avg)
    np.savetxt("%s/Mut_%d/contact_Qi/sim_feature_err.dat" % (name,iteration), sim_feature_err)
    np.savetxt("%s/Mut_%d/contact_Qi/Jacobian.dat" % (name,iteration), Jacobian_avg)
    np.savetxt("%s/Mut_%d/contact_Qi/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)
