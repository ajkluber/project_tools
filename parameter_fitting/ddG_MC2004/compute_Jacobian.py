""" Submodule to calculate the simulation phi values

Description:

    Submodule to calculate the energetic perturbation from each
mutation and the results DeltaDetla G's for simulation.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import os
import sys
import argparse
import time
import numpy as np

from mutatepdbs import get_core_mutations, get_scanning_mutations, get_exp_ddG

#import project_tools.parameter_fitting.util.util as util
from project_tools.parameter_fitting.util.util import *


global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

global SKIP_INTERACTIONS
SKIP_INTERACTIONS = [1,8,9]

def get_dHk(model,rij,Fij_conts,Fij):
    """ Get perturbed potential energy """
    dHk = np.zeros(rij.shape[0],float)
    for i in range(len(Fij_conts)): 
        cont_idx = Fij_conts[i]
        dHk = dHk - Fij[i]*model.pair_eps[cont_idx]*model.pair_V[cont_idx](rij[:,cont_idx])
    return dHk

def get_Vp_plus_Vpk(model,Vp,rij,Fij_conts,Fij):
    """ Get perturbed potential energy """
    Vp_plus_Vpk = np.array(Vp,copy=True)
    for i in range(len(Fij_conts)):
        cont_idx = Fij_conts[i]
        param_idx = model.pairwise_param_assignment[cont_idx]
        Vp_plus_Vpk[:,param_idx] = Vp_plus_Vpk[:,param_idx] - Fij[i]*model.pair_V[cont_idx](rij[:,cont_idx])
    return Vp_plus_Vpk

def get_dHk_for_state(model,rij,Fij_pairs,Fij,state,n_frames):
    """ Get perturbed potential energy """
    dHk_state = np.zeros(n_frames,float)
    for i in range(len(Fij_pairs)):
        pair = Fij_pairs[i]
        
        # Loop over interactions that a given pair have.
        flag = (model.pairs[:,0] == pair[0]).astype(int)*(model.pairs[:,1] == pair[1]).astype(int)
        pair_interaction_indices = np.where(flag == 1)[0]
        for j in range(len(pair_interaction_indices)):
            inter_idx = pair_interaction_indices[j]
            param_idx = model.pairwise_param_assignment[inter_idx]
            # If that interaction corresponds to a fitting parameter 
            # and is not an excluded volume interaction (e.g. LJ12)
            # then add the perturbed interaction energy.
            #if (not (model.pair_type[inter_idx] in SKIP_INTERACTIONS)) and (param_idx in model.fitting_params):
            if param_idx in model.fitting_params:
                dHk_state = dHk_state - Fij[i]*model.pair_eps[inter_idx]*model.pair_V[inter_idx](rij[state,inter_idx])
    return dHk_state

def get_Vp_plus_Vpk_for_state(model,Vp,rij,Fij_pairs,Fij,state):
    """ Get perturbed potential energy """
    Vp_plus_Vpk_state = np.array(Vp,copy=True)
    for i in range(len(Fij_pairs)):     
        pair = Fij_pairs[i]

        # Loop over interactions that a given pair have.
        flag = (model.pairs[:,0] == pair[0]).astype(int)*(model.pairs[:,1] == pair[1]).astype(int)
        pair_interaction_indices = np.where(flag == 1)[0]
        for j in range(len(pair_interaction_indices)):
            inter_idx = pair_interaction_indices[j]
            param_idx = model.pairwise_param_assignment[inter_idx]
            # If that interaction corresponds to a fitting parameter 
            # and is not an excluded volume interaction (e.g. LJ12)
            # then add the perturbed interaction energy.
            #if (not (model.pair_type[inter_idx] in SKIP_INTERACTIONS)) and (param_idx in model.fitting_params):
            if param_idx in model.fitting_params:
                # If that interaction is associated with one of the parameters being fit.
                fitting_param_idx = np.where(model.fitting_params == param_idx)[0][0]
                change = Fij[i]*model.pair_V[inter_idx](rij[state,inter_idx])
                Vp_plus_Vpk_state[:,fitting_param_idx] = Vp_plus_Vpk_state[:,fitting_param_idx] - change

    return Vp_plus_Vpk_state

def get_target_feature(model,fitopts):
    """ Get target features """
    name = model.name
    iteration = fitopts['iteration']
    
    cwd = os.getcwd()
    os.chdir("%s/mutants" % name)
    target_feature, target_feature_err = get_exp_ddG()
    os.chdir(cwd)

    return target_feature, target_feature_err

def calculate_average_Jacobian(model,fitopts,scanning_only=False,scanfij=0.5,saveas="Q_phi.dat",test=False):
    """ Calculate ddG's and Jacobian over replicas"""
    
    name = model.name
    iteration = fitopts['iteration']

    # Get mutations and fraction of native pairs deleted for each mutation.
    os.chdir("%s/mutants" % name)
    mutants = get_core_mutations()
    n_muts = len(mutants)
    Fij, Fij_pairs = get_mutant_fij(model,fitopts,mutants)
    os.chdir("../iteration_%d" % iteration)

    Tlist = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(Tlist[0].split("_")[0]))
    
    # Get pairwise distances from trajectories
    trajfiles = [ "%s/traj.xtc" % x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 
    native = "%s/Native.pdb" % Tlist[0]
    rij = get_rij(model,trajfiles,native)

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_pairs]

    # Found speedups by calculating quantities on per state basis
    U,TS,N,Uframes,TSframes,Nframes = concatenate_state_indicators(Tlist,bounds,coord="Q.dat")

    # Average dimensionless potential energy for each state
    Vp_U   = get_Vp_for_state(model,rij,U,Uframes)
    Vp_TS  = get_Vp_for_state(model,rij,TS,TSframes)
    Vp_N   = get_Vp_for_state(model,rij,N,Nframes)
    sumVp_U = np.mean(Vp_U,axis=0); sumVp_TS = np.mean(Vp_TS,axis=0); sumVp_N = np.mean(Vp_N,axis=0) 

    # Free energy change upon mutation
    dG = np.zeros((3,n_muts),float); ddG = np.zeros((2,n_muts),float); phi = np.zeros(n_muts,float) 

    # Initialize Jacobian. 
    Jacobian = np.zeros((2*n_muts,model.n_fitting_params),float)
    sim_feature = np.zeros(2*n_muts,float)

    lasttime = time.time()
    for k in range(n_muts):
        sys.stdout.write("\r\x1b[K" + " Computing Jacobian: %3f %% done " % (float(k)/float(n_muts)))
        sys.stdout.flush()
        compute_mutation(k,beta,model,rij,n_muts,sim_feature,Jacobian,dG,ddG,phi,Fij_pairs,Fij,
                        U,TS,N,Uframes,TSframes,Nframes,Vp_U,Vp_TS,Vp_N,sumVp_U,sumVp_TS,sumVp_N)

    print "\n  Avg time per mutation: %.2f sec" % ((time.time() - lasttime)/float(n_muts))

    savedir = "newton"
    if test:
        savedir = "test"

    if not os.path.exists("%s" % savedir):
        os.mkdir("%s" % savedir)
    np.savetxt("%s/Jacobian.dat" % savedir,Jacobian)
    np.savetxt("%s/sim_feature.dat" % savedir,sim_feature)
    phi_string = save_phi_values(mutants,state_labels,dG,ddG,phi)
    open("%s/%s" % (savedir,saveas),"w").write(phi_string)

    os.chdir("../..")

    # Call function to compute error on sim_feature.
    Jacobian_err = np.zeros(Jacobian.shape,float)
    sim_feature_err = None

    return sim_feature, sim_feature_err, Jacobian, Jacobian_err

def compute_mutation(k,beta,model,rij,n_muts,sim_feature,Jacobian,dG,ddG,phi,Fij_pairs,Fij,
                U,TS,N,Uframes,TSframes,Nframes,Vp_U,Vp_TS,Vp_N,sumVp_U,sumVp_TS,sumVp_N):

    # Compute energy perturbation
    dHk_U  = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],U,Uframes)
    dHk_TS = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],TS,TSframes)
    dHk_N  = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],N,Nframes)

    # Free energy perturbation formula. Equation (4) in reference (1).
    expdHk_U  = np.mean(np.exp(-beta*dHk_U))
    expdHk_TS = np.mean(np.exp(-beta*dHk_TS))
    expdHk_N  = np.mean(np.exp(-beta*dHk_N))
    dG_U  = -np.log(expdHk_U); dG_TS = -np.log(expdHk_TS); dG_N  = -np.log(expdHk_N)

    # DeltaDeltaG's. Equations (5) in reference (1).
    ddG_stab = (dG_N - dG_U); ddG_dagg = (dG_TS - dG_U)

    dG[0,k] = dG_U; dG[1,k] = dG_TS; dG[2,k] = dG_N
    ddG[0,k] = ddG_dagg; ddG[1,k] = ddG_stab

    #print dG_U, dG_TS, dG_N, ddG_stab, ddG_dagg     # DEBUGGING

    # Phi-value
    if ddG_stab != 0:
        phi[k] = ddG_dagg/ddG_stab
    else:
        phi[k] = 0 

    # simulation feature 
    sim_feature[k] = ddG_dagg
    sim_feature[k + n_muts] = ddG_stab

    # Get perturbed dimensionless potential energy.
    Vp_plus_Vpk_U  = get_Vp_plus_Vpk_for_state(model,Vp_U,rij,Fij_pairs[k],Fij[k],U)
    Vp_plus_Vpk_TS = get_Vp_plus_Vpk_for_state(model,Vp_TS,rij,Fij_pairs[k],Fij[k],TS)
    Vp_plus_Vpk_N  = get_Vp_plus_Vpk_for_state(model,Vp_N,rij,Fij_pairs[k],Fij[k],N)

    # Thermal averages for matrix equation (9).
    # Vectorized computation of Jacobian. 5x faster than a for loop.
    Vp_Vpk_expdHk_U  = sum((Vp_plus_Vpk_U.T*np.exp(-beta*dHk_U)).T)/float(Uframes)
    Vp_Vpk_expdHk_TS = sum((Vp_plus_Vpk_TS.T*np.exp(-beta*dHk_TS)).T)/float(TSframes)
    Vp_Vpk_expdHk_N  = sum((Vp_plus_Vpk_N.T*np.exp(-beta*dHk_N)).T)/float(Nframes)

    Jacobian[k,:] = beta*(((Vp_Vpk_expdHk_TS/expdHk_TS) - sumVp_TS) - ((Vp_Vpk_expdHk_U/expdHk_U) - sumVp_U))
    Jacobian[k + n_muts,:] = beta*(((Vp_Vpk_expdHk_N/expdHk_N) - sumVp_N) - ((Vp_Vpk_expdHk_U/expdHk_U) - sumVp_U))

def compute_Jacobian_for_directory(model,beta,mutants,Fij,Fij_pairs,bounds,state_labels,saveas="Q_phi.dat",test=False):
    """ Calculates the feature vector (ddG's) and Jacobian for one directory """
    # Get trajectory, state indicators, contact energy
    traj,rij = get_traj_rij(model)

    Q = np.loadtxt("Q.dat")
    U,TS,N,Uframes,TSframes,Nframes = get_state_indicators(Q,bounds)

    Vp_U   = get_Vp_for_state(model,rij,U,Uframes)
    Vp_TS  = get_Vp_for_state(model,rij,TS,TSframes)
    Vp_N   = get_Vp_for_state(model,rij,N,Nframes)

    # Average dimensionless potential energy for each state.
    sumVp_U   = np.mean(Vp_U,axis=0) 
    sumVp_TS  = np.mean(Vp_TS,axis=0) 
    sumVp_N   = np.mean(Vp_N,axis=0) 

    # Compute deltaG for each state. Then DeltaDelta G with respect to the
    # first state (assumed to be the unfolded/denatured state).
    # Units of kT.
    dG = np.zeros((3,len(mutants)),float)
    ddG = np.zeros((2,len(mutants)),float)
    phi = np.zeros(len(mutants),float) 

    # Initialize Jacobian
    #Jacobian = np.zeros((2*len(mutants),model.n_model_param),float)
    Jacobian = np.zeros((2*len(mutants),model.n_fitting_params),float)
    sim_feature = np.zeros(2*len(mutants),float)

    avg_rowtime = np.zeros(len(mutants),float)
    lasttime = time.time()
    for k in range(len(mutants)):
        mut = mutants[k]
        # Compute energy perturbation
        dHk_U  = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],U,Uframes)
        dHk_TS = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],TS,TSframes)
        dHk_N  = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],N,Nframes)

        # Free energy perturbation formula. Equation (4) in reference (1).
        expdHk_U  = np.mean(np.exp(-beta*dHk_U))
        expdHk_TS = np.mean(np.exp(-beta*dHk_TS))
        expdHk_N  = np.mean(np.exp(-beta*dHk_N))
        dG_U  = -np.log(expdHk_U)
        dG_TS = -np.log(expdHk_TS)
        dG_N  = -np.log(expdHk_N)

        # DeltaDeltaG's. Equations (5) in reference (1).
        ddG_stab = (dG_N - dG_U)
        ddG_dagg = (dG_TS - dG_U)

        #print dG_U, dG_TS, dG_N, ddG_stab, ddG_dagg     # DEBUGGING

        # Phi-value
        if ddG_stab != 0:
            phi_value = ddG_dagg/ddG_stab
        else:
            phi_value = 0

        dG[0,k] = dG_U
        dG[1,k] = dG_TS
        dG[2,k] = dG_N
        ddG[0,k] = ddG_dagg
        ddG[1,k] = ddG_stab
        phi[k] = phi_value

        # simulation feature 
        sim_feature[k] = ddG_dagg
        sim_feature[k + len(mutants)] = ddG_stab

        # Get perturbed dimensionless potential energy.
        Vp_plus_Vpk_U  = get_Vp_plus_Vpk_for_state(model,Vp_U,rij,Fij_pairs[k],Fij[k],U)
        Vp_plus_Vpk_TS = get_Vp_plus_Vpk_for_state(model,Vp_TS,rij,Fij_pairs[k],Fij[k],TS)
        Vp_plus_Vpk_N  = get_Vp_plus_Vpk_for_state(model,Vp_N,rij,Fij_pairs[k],Fij[k],N)

        # Thermal averages for matrix equation (9).
        # Vectorized computation of Jacobian. 5x faster than a for loop.
        Vp_Vpk_expdHk_U  = sum((Vp_plus_Vpk_U.T*np.exp(-beta*dHk_U)).T)/Uframes
        Vp_Vpk_expdHk_TS = sum((Vp_plus_Vpk_TS.T*np.exp(-beta*dHk_TS)).T)/TSframes
        Vp_Vpk_expdHk_N  = sum((Vp_plus_Vpk_N.T*np.exp(-beta*dHk_N)).T)/Nframes

        Jacobian[k,:] = beta*(((Vp_Vpk_expdHk_TS/expdHk_TS) - sumVp_TS) - ((Vp_Vpk_expdHk_U/expdHk_U) - sumVp_U))
        Jacobian[k + len(mutants),:] = beta*(((Vp_Vpk_expdHk_N/expdHk_N) - sumVp_N) - ((Vp_Vpk_expdHk_U/expdHk_U) - sumVp_U))

        thistime = time.time()
        dt = thistime - lasttime
        avg_rowtime[k] = dt
        lasttime = thistime
        print "    mutant %d  %5s    %.2f seconds = %.2f min" % (k,mut,dt,dt/60.)

    print "  Avg rowtime: %.2f sec" % np.mean(avg_rowtime)
    if test:
        if not os.path.exists("test"):
            os.mkdir("test")
        np.savetxt("test/Jacobian.dat",Jacobian)
        np.savetxt("test/sim_feature.dat",sim_feature)
        phi_string = save_phi_values(mutants,state_labels,dG,ddG,phi)
        open("test/%s" % saveas,"w").write(phi_string)
    else:
        if not os.path.exists("mut"):
            os.mkdir("mut")
        np.savetxt("mut/Jacobian.dat",Jacobian)
        np.savetxt("mut/sim_feature.dat",sim_feature)
        phi_string = save_phi_values(mutants,state_labels,dG,ddG,phi)
        open("mut/%s" % saveas,"w").write(phi_string)

    return sim_feature, Jacobian

def compute_dHk(model,fitopts):

    name = model.name
    iteration = fitopts['iteration']

    # Get list of mutations and fraction of native pairs deleted for 
    # each mutation.
    os.chdir("%s/mutants" % name)
    mutants = get_core_mutations()
    Fij, Fij_pairs = get_mutant_fij(model,fitopts,mutants)
    os.chdir("../..")

    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)

    dirs = [ x.rstrip("\n") for x in open("long_temps","r").readlines() ]
    temps = [ x.split('_')[0] for x in dirs ]

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_pairs]

    # Loop over temperatures.
    for n in range(len(dirs)):
        T = temps[n]
        dir = dirs[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        os.chdir(dir)
        if not os.path.exists("expdHk_%s_U.dat" % mutants[0]):
            print "going into directory: %s" % dir
            # Get trajectory, state indicators, contact energy
            traj,rij = get_traj_rij(model)

            Q = np.loadtxt("Q.dat")
            U,TS,N,Uframes,TSframes,Nframes = get_state_indicators(Q,bounds)
            #allframes = np.ones(len(traj.n_frames)).astype(bool)

            for k in range(len(mutants)):
            #for k in range(len(mutants)):
                mut = mutants[k]
                # Compute energy perturbation
                #dHk = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],allframes,traj.n_frames)
                #np.savetxt("dHk_%s_U.dat" % mut,dHk_U)

                dHk_U = np.zeros(traj.n_frames,float)
                dHk_TS = np.zeros(traj.n_frames,float)
                dHk_N = np.zeros(traj.n_frames,float)
                dHk_U[U] = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],U,Uframes)
                dHk_TS[TS] = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],TS,TSframes)
                dHk_N[N] = get_dHk_for_state(model,rij,Fij_pairs[k],Fij[k],N,Nframes)

                np.savetxt("dHk_%s_U.dat" % mut,dHk_U)
                np.savetxt("dHk_%s_TS.dat" % mut,dHk_TS)
                np.savetxt("dHk_%s_N.dat" % mut,dHk_N)

                print "  mutant %d of %d done %s " % (k,len(mutants),mut)

        os.chdir("..")

    os.chdir("../..")

def get_mutant_fij(model,fitopts,mutants):
    """ Load in the fraction of contact loss for each mutation.

    Description:

        Since the fraction of contacts lost matrix f^k_ij is sparse only load
    in the of nonzero elements and their indices. Determine the contact indices
    for nonzero entries. Let only allow fij's for contacts. 
    """
    Fij_pairs = []
    Fij = []
    for mut in mutants:
        if fitopts["nonnative"]:
            fij_temp = np.loadtxt("fij_%s.dat" % mut)
        else:
            fij_temp = model.Qref*np.loadtxt("fij_%s.dat" % mut)
        indices = np.nonzero(fij_temp)
        Fij.append(fij_temp[indices])
        #mean_fij = np.mean(fij_temp[indices])
        #Fij.append(mean_fij*np.ones(len(fij_temp[indices])))        # DEBUGGING
        Fij_pairs.append(np.array(zip(indices[0]+1,indices[1]+1)))
        
    return Fij, Fij_pairs

def get_mutant_fij_scanning(model, mutants, fij=0.5):
    """Estimate the local contact fraction loss for scanning mutations.
    
        The ddGs for scanning mutations could be affected by a multiplicity of
    factors that exceed the mere contacts lost (e.g. effect of solvent,
    interaction with helix dipole, interaction with charged residues, loss of
    hydrogen bonds, among others). As a result, the exact determination varies
    according to which factor(s) weigh the most. Given that this is a
    coarse-grained model, as a first approach we will estimate an average value
    of fij = 0.5. This is as said before arbitrary, but should give us an
    estimate of how the local epsilons vary according to the value of the ddGs
    involved.

    Also, on a first approach only the [i, i+4] contacts will be affected by
    the said value of fij.

    """
    
    Fij_pairs = []
    Fij_conts = []
    Fij = []

    if mutants == []:
        pass
    else:
        for mut in mutants:
            # For the time being only keeping i, i+4
            mut_res_number = int("".join(list(mut)[1:-1]))
            Fij.append(fij)
            temppairs = []
            tempconts = []
            tempfij = []

            for j in range(model.n_pairs):
                if (model.pairs[j,0] == mut_res_number) and (model.pairs[j,1] == (mut_res_number+4)):
                    contact_num = j
                    temppairs.append([mut_res_number-1,mut_res_number+3])  #i, i+4
                    tempconts.append(contact_num)
                    tempfij.append(fij)
                    break
                else:
                    continue

            Fij_conts.append(np.array(tempconts))
            Fij_pairs.append(temppairs)
    return Fij, Fij_pairs, Fij_conts

def save_phi_values(mutants,state_labels,dG,ddG,phi):
    """ Save the calculated dG, ddG, and phi values for states"""

    header_string = "# mut" 
    for i in range(len(state_labels)):
        header_string += "     dG_%s(kT)" % state_labels[i]
    header_string += 4*" "
    for i in range(1,len(state_labels)):
        header_string += " ddG_%s-%s(kT)" % (state_labels[i],state_labels[0])
    for i in range(1,len(state_labels)-1):
        header_string += "   Phi_%s/%s" % (state_labels[i],state_labels[-1])

    data_string = ''
    for j in range(len(mutants)):
        line = "%6s"%mutants[j]
        for i in range(len(dG)):
            line += "  %10.5f " % dG[i][j]
        for i in range(len(ddG)):
            line += "  %10.5f " % ddG[i][j]
        line += "  %10.5f  " % phi[j]
        data_string += line+"\n"
    print "ddG and Phi values:"
    print header_string
    print data_string
    return header_string+"\n"+data_string

if __name__ == "__main__":
    import model_builder as mdb

    #parser = argparse.ArgumentParser(description='Calculate .')
    #parser.add_argument('--name', type=str, required=True, help='name.')
    #parser.add_argument('--iteration', type=int, required=True, help='iteration.')
    #args = parser.parse_args()
    
    #name = args.name
    #iteration= args.iteration
    
    name = "S6"
    model, fitopts = mdb.inputs.load_model(name)
    #compute_dHk(model,fitopts)
    sim_f, sim_f_err, J, J_err = calculate_average_Jacobian(model,fitopts)
    target_feature, target_feature_err = get_target_feature(model,fitopts)
    np.savetxt("S6/iteration_0/newton/target_feature.dat", target_feature)
