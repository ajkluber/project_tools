

import os
import argparse
import time
import numpy as np

import mdtraj as md

from mutatepdbs import get_core_mutations, get_scanning_mutations, get_exp_ddG
from compute_Jacobian import get_dHk_for_state,get_mutant_fij,save_phi_values

#import project_tools.parameter_fitting.util.util as util
from project_tools.parameter_fitting.util.util import *


global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621


def calculate_perturbation_phi(model,fij):
    ''' Calculate the average feature vector (ddG's) and Jacobian '''
    
    name = model.subdir
    iteration = model.iteration

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    os.chdir("%s/mutants" % name)
    ## Get list of mutations and fraction of native contacts deleted for 
    ## each mutation.
    mutants = get_core_mutations()

    Fij, Fij_pairs, Fij_conts = get_mutant_fij(model,mutants)
    Fij = [ fij*np.ones(len(x)) for x in Fij ]

    os.chdir(sub)
    temperatures = [ x.split('_')[0] for x in open("long_temps_last","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in iteration subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. Save the average
    ## feature vector and Jacobian in the iteration/newton directory.
    phi_all = []
    ddG_all = []
    lasttime = time.time()
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for iteration_%d/%s" % (model.iteration,dir)
        os.chdir(dir)
        phi,ddG = compute_perturbation_phi_for_directory(model,beta,mutants,Fij,Fij_pairs,Fij_conts,bounds,state_labels)
        phi_all.append(phi)
        ddG_all.append(ddG)
        os.chdir("..")
        thistime = time.time()
        timediff = thistime - lasttime
        lasttime = thistime
        print "  calculation took %.2f seconds = %.2f minutes" % (timediff,timediff/60.)

    phi_all = np.array(phi_all)
    ddG_all = np.array(ddG_all)

    ## Take avg. and use standard deviation as error bars.
    phi_avg = sum(phi_all)/float(len(directories))
    phi_err = np.std(phi_all,axis=0)
    ddG_avg = sum(ddG_all)/float(len(directories))
    ddG_err = np.std(ddG_all,axis=0)


    if not os.path.exists("test"):
        os.mkdir("test")
    np.savetxt("test/phi",np.array(zip(phi_avg,phi_err)))
    np.savetxt("test/ddG",np.array(zip(ddG_avg,ddG_err)))

    os.chdir(cwd)

    return phi_avg, phi_err, ddG_avg, ddG_err

def compute_perturbation_phi_for_directory(model,beta,mutants,Fij,Fij_pairs,Fij_conts,bounds,state_labels):
    ''' Calculates the feature vector (ddG's) and Jacobian for one directory '''
    ## Get trajectory, state indicators, contact energy
    traj,rij = get_traj_rij(model)

    Q = np.loadtxt("Q.dat")
    U,TS,N,Uframes,TSframes,Nframes = get_state_indicators(Q,bounds)

    ## Compute deltaG for each state. Then DeltaDelta G with respect to the
    ## first state (assumed to be the unfolded/denatured state).
    ## Units of kT.
    dG = np.zeros((3,len(mutants)),float)
    #ddG = np.zeros((2,len(mutants)),float)
    ddG = np.zeros(2*len(mutants),float)
    phi = np.zeros(len(mutants),float) 

    avg_rowtime = np.zeros(len(mutants),float)
    lasttime = time.time()
    for k in range(len(mutants)):
        mut = mutants[k]
        ## Compute energy perturbation
        dHk_U  = get_dHk_for_state(model,rij,Fij_conts[k],Fij[k],U,Uframes)
        dHk_TS = get_dHk_for_state(model,rij,Fij_conts[k],Fij[k],TS,TSframes)
        dHk_N  = get_dHk_for_state(model,rij,Fij_conts[k],Fij[k],N,Nframes)

        ## Free energy perturbation formula. Equation (4) in reference (1).
        expdHk_U  = np.mean(np.exp(-beta*dHk_U))
        expdHk_TS = np.mean(np.exp(-beta*dHk_TS))
        expdHk_N  = np.mean(np.exp(-beta*dHk_N))
        dG_U  = -np.log(expdHk_U)
        dG_TS = -np.log(expdHk_TS)
        dG_N  = -np.log(expdHk_N)

        ## DeltaDeltaG's. Equations (5) in reference (1).
        ddG_stab = (dG_N - dG_U)
        ddG_dagg = (dG_TS - dG_U)

        ## Phi-value
        phi_value = ddG_dagg/ddG_stab

        ddG[k] = ddG_dagg
        ddG[k+len(mutants)] = ddG_stab
        phi[k] = phi_value

        thistime = time.time()
        dt = thistime - lasttime
        avg_rowtime[k] = dt
        lasttime = thistime
        print "    mutant %d  %5s    %.2f seconds = %.2f min" % (k,mut,dt,dt/60.)

    print "  Avg rowtime: %.2f sec" % np.mean(avg_rowtime)

    return phi,ddG


if __name__ == "__main__":
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='name.')
    parser.add_argument('--iteration', type=int, required=True, help='iteration.')
    args = parser.parse_args()
    
    name = args.name
    iteration= args.iteration

    fij = 0.5

    contacts = np.loadtxt("%s/contacts.dat" % name)
    pdb = "%s.pdb" % name
    defaults = True
    model = mdb.models.SmogCalpha.SmogCalpha(pdb=pdb,contacts=contacts,defaults=defaults,iteration=iteration)
    phi_avg, phi_err, ddG_avg, ddG_err = calculate_perturbation_phi(model,fij)

