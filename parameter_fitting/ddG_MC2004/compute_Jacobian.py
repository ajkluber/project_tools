""" Submodule to calculate the simulation phi values

Description:

    Submodule to calculate the energetic perturbation from each
mutation and the results DeltaDetla G's for simulation.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
"""

import numpy as np
import os
import argparse

import mdtraj as md

import model_builder.models as models

from mutatepdbs import get_core_mutations, get_exp_ddG

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_average_Jacobian(model,append_log):
    """ Calculate the average feature vector (ddG's) and Jacobian """
    
    append_log(model.subdir,"Starting: Calculating_Jacobian")

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)
    os.chdir(model.subdir+"/mutants")
    ## Get list of mutations and fraction of native contacts deleted for 
    ## each mutation.
    mutants = get_core_mutations()
    ddGexp, ddGexp_err = get_exp_ddG()
    Fij, Fij_pairs, Fij_conts = get_mutant_fij(model,mutants)
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    epsilons = model.contact_epsilons
    deltas = model.contact_deltas
    sigmas = model.contact_sigmas

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in Mut subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. Save the average
    ## feature vector and Jacobian in the Mut/newton directory.
    sim_feature_all = []
    Jacobian_all = []
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for Mut_%d/%s" % (model.Mut_iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,traj,dir,beta,mutants,Fij_pairs,Fij_conts,bounds,state_labels,epsilons,delta,sigmas)
        sim_feature_all.append(sim_feature)
        Jacobian_all.append(Jacobian)
        os.chdir("..")

    sim_feaure_all = np.array(sim_feaure_all)
    Jacobian_all = np.array(Jacobian_all)

    ## Take avg. and use standard deviation as error bars.
    sim_feature_avg = sum(sim_feature_all)/float(len(directories))
    sim_feature_err = np.std(sim_feature_all,axis=0)
    Jacobian_avg = sum(Jacobian_all)/float(len(directories))
    Jacobian_err = np.std(Jacobian_all,axis=0)

    if not os.path.exists("newton"):
        os.mkdir("newton")

    print "  Saving feature vector and Jacobian in Mut_%d/newton" % model.Mut_iteration
    np.savetxt("newton/target_feature.dat",ddGexp)
    np.savetxt("newton/target_feature_err.dat",ddGexp_err)
    np.savetxt("newton/sim_feature.dat",sim_feature_avg)
    np.savetxt("newton/sim_feature_err.dat",sim_feature_err)
    np.savetxt("newton/Jacobian.dat",Jacobian_avg)
    np.savetxt("newton/Jacobian_err.dat",Jacobian_err)

    os.chdir(cwd)
    append_log(model.subdir,"Finished: Calculating_Jacobian")

def compute_Jacobian_for_directory(model,traj,dir,beta,mutants,Fij_pairs,Fij_conts,bounds,state_labels,epsilons,delta,sigmas):
    """ Calculates the feature vector (ddG's) and Jacobian for one directory """
    ## Get trajectory, state indicators, contact energy
    traj,U,TS,N,Uframes,TSframes,Nframes,Vij = get_states_Vij(model,bounds,epsilons,delta,sigmas)

    ## Avg. contact energy for each state.
    Vij_U  = sum(Vij[U,:])/Uframes
    Vij_TS = sum(Vij[TS,:])/TSframes
    Vij_N  = sum(Vij[N,:])/Nframes

    ## Compute deltaG for each state. Then DeltaDelta G with respect to the
    ## first state (assumed to be the unfolded/denatured state).
    ## Units of kT.
    dG = [[],[],[]]
    ddG = [[],[]]
    phi = []

    ## Initialize Jacobian
    Jacobian = np.zeros((2*len(mutants),model.n_contacts),float)
    sim_feature = np.zeros(2*len(mutants),float)
    

    for k in range(len(mutants)):
        mut = mutants[k]
        print "    row %d   mutant %s" % (k,mut)
        ## Compute energy perturbation
        dHk = compute_dHk(model,traj,mut,Fij_pairs[k],Fij_conts[k])

        ## Free energy perturbation formula. Equation (4) in reference (1).
        dG_U  = -np.log(np.sum(np.exp(-beta*dHk[U]))/Uframes)
        dG_TS = -np.log(np.sum(np.exp(-beta*dHk[TS]))/TSframes)
        dG_N  = -np.log(np.sum(np.exp(-beta*dHk[N]))/Nframes)

        ## DeltaDeltaG's. Equations (5) in reference (1).
        ddG_stab = (dG_N - dG_U)
        ddG_dagg = (dG_TS - dG_U)

        ## Phi-value
        phi_value = ddG_dagg/ddG_stab

        dG[0].append(dG_U)
        dG[1].append(dG_TS)
        dG[2].append(dG_N)
        ddG[0].append(ddG_dagg)
        ddG[1].append(ddG_stab)
        phi.append(phi_value)

        ## Thermal averages for matrix equation (9).
        expdHk_U  = sum(np.exp(-beta*dHk[U]))/Uframes
        expdHk_TS = sum(np.exp(-beta*dHk[TS]))/TSframes
        expdHk_N  = sum(np.exp(-beta*dHk[N]))/Nframes

        Vij_expdHk_U  = sum((Vij[U,:].T*np.exp(-beta*dHk[U])).T)/Uframes
        Vij_expdHk_TS = sum((Vij[TS,:].T*np.exp(-beta*dHk[TS])).T)/TSframes
        Vij_expdHk_N  = sum((Vij[N,:].T*np.exp(-beta*dHk[N])).T)/Nframes

        ## Compute all columns with Fij_k zero.
        Jacobian[k,:] = -beta*((Vij_TS - Vij_U) -  ((Vij_expdHk_TS/expdHk_TS) - (Vij_expdHk_U/expdHk_U)))
        Jacobian[k + len(mutants),:] = -beta*((Vij_N - Vij_U)  -  ((Vij_expdHk_N/expdHk_N) - (Vij_expdHk_U/expdHk_U)))

        ## Replace columns for which Fij_k is not zero.
        Jacobian[k,Fij_conts[k]] = -beta*((Vij_TS[Fij_conts[k]] - Vij_U[Fij_conts[k]])  - \
            (1. - Fij[k])*((Vij_expdHk_TS[Fij_conts[k]]/expdHk_TS) - \
                           (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))

        Jacobian[k + len(mutants),Fij_conts[k]] = -beta*((Vij_N[Fij_conts[k]] - Vij_U[Fij_conts[k]])  -  \
            (1. - Fij[k])*((Vij_expdHk_N[Fij_conts[k]]/expdHk_N)   - \
                          (Vij_expdHk_U[Fij_conts[k]]/expdHk_U)))

        ## return feature vector and Jacobian
        sim_feature[:len(mutants)] = ddG_dagg
        sim_feature[len(mutants):] = ddG_stab

    np.savetxt("mut/Jacobian.dat",Jacobian)
    np.savetxt("mut/sim_feature.dat",sim_feature)
    save_phi_values(mutants,"Q",state_labels,dG,ddG,phi)

    return sim_feature, Jacobian

def get_states_Vij(model,bounds,epsilons,delta,sigmas):
    """ Load trajectory, state indicators, and contact energy """

    traj = md.load("traj.xtc",top="Native.pdb")     ## Loading from file takes most time.
    rij = md.compute_distances(traj,model.contacts-np.ones(model.contacts.shape))
    Q = np.loadtxt("Q.dat")

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

    ## Only count values of potential energy function where interaction is
    ## attractive.
    x = sigmas/rij
    x[(x > 1.09)] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero. 
    Vij = epsilons*(5.*(x**12) - 6.*deltas*(x**10))     ## To Do: Generalize to other contact functions

    return traj,U,TS,N,Uframes,TSframes,Nframes,Vij

def compute_dHk(model,traj,mut,pairs,conts):
    """ Calculate energetic perturbation to each frame from mutation """

    ## Use mdtraj to compute the distances between pairs.
    rij = md.compute_distances(traj,pairs)
    
    ## Epsilons, deltas, and sigmas for relevant pairs for this mutation.
    eps = model.contact_epsilons[conts]
    deltas = model.contact_deltas[conts]
    sigmas = model.contact_sigmas[conts]

    ## Sometimes the energies calculated here blow up because small
    ## changes in distance can lead to large energies on the
    ## repulsive wall. (The actual energies of the simulation do not
    ## blow up). Since this calculation is meant to estimation of 
    ## the change in energy upon mutation, dH_k, I just use frames
    ## where the contact energy LJ12-10 is negative. Without this
    ## consideration, the results are crazy
    x = sigmas/rij
    x[(x > 1.09)] = 1.09  # <-- 1.09 is where LJ12-10 crosses zero.
    #x[(x > 1.211)] = 1.211     ## Old, arbitrary cutoff.

    ## Calculate dH_k using distances and parameters. Save.
    Vij = -Fij[k]*eps*(5.*(x**12) - 6.*deltas*(x**10))
    dH_k = sum(Vij.T)
    np.savetxt("dH_"+mut+".dat",dH_k)

    return dH_k

def get_mutant_fij(model,mutants):
    """ Load in the fraction of contact loss for each mutation.

    Description:

        Since the fraction of contacts lost matrix f^k_ij is sparse only load
    in the of nonzero elements and their indices. Determine the contact indices
    for nonzero entries.
    """
    Fij_pairs = []
    Fij_conts = []
    Fij = []
    for mut in mutants:
        fij_temp = model.Qref*np.loadtxt("fij_"+mut+".dat")
        indices = np.nonzero(fij_temp)
        Fij.append(fij_temp[indices])
        temppairs = []
        tempconts = []
        tempfij = []
        for i in range(len(indices[0])):
            for j in range(model.n_contacts):
                if (model.contacts[j,0] == (indices[0][i]+1)) and (model.contacts[j,1] == (indices[1][i]+1)):
                    contact_num = j
                    temppairs.append([indices[0][i],indices[1][i]])
                    tempconts.append(contact_num)
                    tempfij.append(fij_temp[indices[0][i],indices[1][i]])
                    break
                else:
                    continue
        Fij_conts.append(np.array(tempconts))
        Fij_pairs.append(temppairs)
    return Fij, Fij_pairs, Fij_conts

def get_Qij(model,r,sig,delta,interaction_nums):
    """ Calculates the normalized interaction betwen nonbonded pairs.

        Might use to generalize for different types of contacts.
    """
    print "  Calculating Qij..."
    qij = model.nonbond_interaction(r,sig,delta)
    return qij

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
        state_bounds.append(float(info[1]))
        state_bounds.append(float(info[2]))
        state_labels.append(info[0])
    
    return state_bounds,state_labels

def save_phi_values(mutants,coord,state_labels,dG,ddG,phi):
    """ Save the calculated dG, ddG, and phi values for states"""

    header_string = "# mut" 
    for i in range(len(state_labels)):
        header_string += "     dG_"+state_labels[i]+"(kT)"
    header_string += "    "
    for i in range(1,len(state_labels)):
        header_string += " ddG_"+state_labels[i]+"-"+state_labels[0]+"(kT)"
    for i in range(1,len(state_labels)-1):
        header_string += "   Phi_"+state_labels[i]+"/"+state_labels[-1]

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

    outputfile = open("phi/"+coord+"_phi.dat","w")
    outputfile.write(header_string+"\n"+data_string)
    outputfile.close()

if __name__ == "__main__":
    ## To Do: Make into a command line utility

    #parser = argparse.ArgumentParser(description='Calculate .')
    #parser.add_argument('--subdir', type=str, required=True, help='Directory.')
    #parser.add_argument('--calc_dH', action='store_true', help='Calculate dH for mutants.')
    #args = parser.parse_args()
    #pdb = args.subdir+".pdb"
    
    def dummy(this,that):
        pass

    pdb = "r15.pdb"
    model = models.SmogCalpha.SmogCalpha(pdb)
    calculate_average_Jacobian(model,append_log)
