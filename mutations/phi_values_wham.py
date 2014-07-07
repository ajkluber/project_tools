
import numpy as np
import os
import argparse

import mdtraj as md

import model_builder.models as models
import model_builder.systems as systems

from mutatepdbs import get_core_mutations

from project_tools.analysis import wham


global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_phi_values_wham(model,append_log):
    """ Calculate the simulation delta G's, deltadelta G's, and phi-values using WHAM. 

    Description:
        Requires the definition of states along Q and the energetic
    perturbation for each mutation dH_k. Saves output in Mut_#/phi

    NOT DONE 7-10-14
    """
    
    #append_log(model.subdir,"Starting: Calculating_phi_values")
    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)


    #beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
    os.chdir(model.subdir)
    if not os.path.exists(sub+"/phi"):
        os.mkdir(sub+"/phi")

    ## Get the mutations.
    os.chdir("mutants")
    mutants = get_core_mutations()
    os.chdir(sub)

    ## Load:
    ## These need to be created by the user after the user messes
    ## with the output free energy and finds where the two basins
    ## are equal.
    ## output temperature <-- sub/whamQ/Tf.txt   
    ## state boundaries   <-- sub/whamQ/state_bounds.txt
    Tf = get_Tf_choice()
    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## loop over mutants
    for k in range(len(mutants)):
    #for k in [0]:       ## DEBUGGING
        mut = mutants[k]
        print "  Running wham for exp(-beta*dH_"+mut+")"
        ## Run WHAM to get <exp(-beta*dH_k)>_X for each mutant.
        ## 
        wham.run_wham_expdH_k(mut,Tf,bounds)

    print "Success!"
    raise SystemExit
    num_states = len(states)
    print "  Loading dH for mutants"

    ## Compute deltaG for each state. Then DeltaDelta G with 
    ## respect to the first state (assumed to be the denatured state).
    ## Units of kT.
    print "  Computing ddG and phi values..."
    dG = [ -np.log(sum(np.exp(-beta*dH[:,states[X]]).T)/float(len(dH[:,states[X]]))) for X in range(num_states) ]
    ddG = [ dG[X]-dG[0] for X in range(1,num_states) ]

    ## Compute the phi value for each mutation. Phi is the ratio
    ## of DeltaDeltaG of the transition state(s) to DeltaDelta G
    ## of the native state (assumed to be the last state). unitless.
    phi = [ ddG[X]/ddG[-1] for X in range(len(ddG)-1) ]
    
    save_phi_values(savedir,mutants,coord,bounds,dG,ddG,phi)
    os.chdir(cwd)
    append_log(model.subdir,"Finished: Calculating_phi_values")
