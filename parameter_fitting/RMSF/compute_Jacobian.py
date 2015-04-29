""" Compute Jacobian for matching a distance distribution 


Description:

    This module computes the jacobian of a distance distribution
such as measured with FRET.

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

def get_rij_Vij(model):
    """ Load trajectory, state indicators, and contact energy """

    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape),periodic=False)
    Vij = model.calculate_contact_potential(rij)

    return traj,rij,Vij

def get_target_feature(model):
    """ Get target features """
    name = model.name
    iteration = model.iteration
    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    

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

def calculate_average_Jacobian(model,residues):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    
    name = model.name
    iteration = model.iteration

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("long_temps_last","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 

    #bounds, state_labels = get_state_bounds()
    #bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in iteration subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. 
    sim_feature_all = []
    Jacobian_all = []
    lasttime = time.time()
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for iteration_%d/%s" % (model.iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,beta,residues)
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

def compute_Jacobian_for_directory(model,beta,residues):
    """ Calculates the feature vector (ddG's) and Jacobian for one directory """

    ## Get trajectory, state indicators, contact energy
    traj,rij,Vij = get_Vij(model)

    ## Get simulation feature
    FRETr = md.compute_distances(traj,residues,periodic=False)
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
    parser.add_argument('--temp', type=str, required=True, help='Temperature.')
    parser.add_argument('--savein', type=str, required=True, help='Iteration.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    temp = args.temp
    savein = args.savein

    model = mdb.check_inputs.load_model(name) 
    model.iteration = iteration
    model.Tf_iteration = iteration

    R_KJ_MOL = 0.0083145
    T = float(temp.split("_")[0])
    beta = 1./(R_KJ_MOL*T)

    cwd = os.getcwd()
    os.chdir("%s/Tf_%d/%s" % (name,iteration,temp))

    #sim_feature, Jacobian = compute_Jacobian_for_directory(model,beta,residues)

    n_bins = 30

    ## Get trajectory, state indicators, contact energy
    traj,rij,Vij = get_rij_Vij(model)

    ## Get simulation feature
    FRETr = (md.compute_distances(traj,residues))[:,0]
    P_r, bins = np.histogram(FRETr,bins=n_bins,density=True)
    
    ## Indicates what bins fall in what bin
    indicator = []
    binframes = []
    print "computing indicator.." 
    for i in range(len(bins)-1):
        temp = ((FRETr > bins[i]).astype(int)*(FRETr <= bins[i+1]).astype(int)).astype(bool)
        indicator.append(temp)
        binframes.append(float(sum((temp == True).astype(int))))

    #print sum(binframes), len(FRETr)       # DEBUGGING

    ## Initialize Jacobian
    Jacobian = np.zeros((n_bins,model.n_contacts),float)

    avg_rij = sum(rij[:,:])/float(traj.n_frames)
    avg_Vij = sum(Vij[:,:])/float(traj.n_frames)

    ## Compute rows of the Jacobian which are correlation functions of 
    ## contact formation with contact energy.
    print "computing Jacobian..."
    for i in range(n_bins):
        #if (i % 10) == 0:
        print "    row %d out of %d" % (i+1,n_bins)
        avg_Vij_r = sum((Vij[indicator[i],:].T).T)/binframes[i]
        Jacobian[i,:] = -beta*(avg_Vij_r - P_r[i]*avg_Vij)

    os.chdir(cwd)

    #target_feature, target_feature_err = get_target_feature(model)
    #sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model,residues)

    method = args.savein
    if not os.path.exists("%s/Tf_%d/%s" % (name,iteration,method)):
        os.mkdir("%s/Tf_%d/%s" % (name,iteration,method))
    np.savetxt("%s/Tf_%d/%s/Jacobian.dat" % (name,iteration,method), Jacobian)
    np.savetxt("%s/Tf_%d/%s/sim_feature.dat" % (name,iteration,method), P_r)
    np.savetxt("%s/Tf_%d/%s/bins.dat" % (name,iteration,method), bins)

    #np.savetxt("%s/Tf_%d/%s/target_feature.dat" % (name,iteration), target_feature)
    #np.savetxt("%s/Tf_%d/%s/target_feature_err.dat" % (name,iteration), target_feature_err)
    #np.savetxt("%s/Tf_%d/%s/target_feature.dat" % (name,iteration), target_feature)
    #np.savetxt("%s/Tf_%d/%s/target_feature_err.dat" % (name,iteration), target_feature_err)
    #np.savetxt("%s/Tf_%d/%s/sim_feature.dat" % (name,iteration), sim_feature_avg)
    #np.savetxt("%s/Tf_%d/%s/sim_feature_err.dat" % (name,iteration), sim_feature_err)
    #np.savetxt("%s/Tf_%d/%s/Jacobian.dat" % (name,iteration), Jacobian_avg)
    #np.savetxt("%s/Tf_%d/%s/Jacobian_err.dat" % (name,iteration) ,Jacobian_err)
