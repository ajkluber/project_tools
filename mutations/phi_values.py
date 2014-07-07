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
import model_builder.systems as systems

from mutatepdbs import get_core_mutations

from project_tools.analysis import wham


global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_dH_for_mutants(model,append_log):
    """ Calculate Hamiltonian perturbation for each mutation dH_k """
    
    append_log(model.subdir,"Starting: Calculating_dH")

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)
    ## Get the fraction of native contacts deleted for each mutation.
    os.chdir(model.subdir+"/mutants")
    mutants = get_core_mutations()
    Fij, Fij_pairs, Fij_conts = get_mutant_fij(model,mutants)
    os.chdir(sub)

    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    ## Loop over all directories in Mut subdirectory. For each trajectory
    ## and each mutation crunch dH.
    for dir in directories:
        os.chdir(dir)
        print "  Loading traj for: ", "Mut_"+str(model.Mut_iteration)+"/"+dir
        traj = md.load("traj.xtc",top="Native.pdb")
        for k in range(len(mutants)):
            mut = mutants[k]
            #if not os.path.exists("dH_"+mut+".dat"):
            if True:
                print "    calc dH for ", mutants[k]
                ## Use mdtraj to compute the distances between pairs.
                rij = md.compute_distances(traj,Fij_pairs[k])
                
                ## Epsilons, deltas, and sigmas for relevant pairs for this mutation.
                eps = model.contact_epsilons[Fij_conts[k]]
                deltas = model.contact_deltas[Fij_conts[k]]
                sigmas = model.contact_sigmas[Fij_conts[k]]

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
            else:
                pass
        os.chdir("..")
    os.chdir(cwd)
    append_log(model.subdir,"Finished: Calculating_dH")

def calculate_phi_values(model,append_log,coord):
    """ Calculate the simulation delta G's, deltadelta G's, and phi-values. 

    Description:

        Requires the definition of states along Q and the energetic
    perturbation for each mutation dH_k. Saves output in each temperature
    directory under phi/ subdirectory.
    """
    
    append_log(model.subdir,"Starting: Calculating_phi_values")
    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)
    
    os.chdir(model.subdir)

    os.chdir("mutants")
    mutants = get_core_mutations()
    os.chdir(sub)

    temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    #T,Tnum = get_Tf_choice()
    #Tdir = "%.2f_%d" % (T,Tnum)

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]
    for n in range(len(directories)):
        T = temperatures[n]
        Tdir = directories[n]
        os.chdir(Tdir)
        print "  Computing ddG and phi-values for ", Tdir
        calculate_ddG_for_temperature(mutants,T,bounds,state_labels)
        os.chdir("..")

    os.chdir(cwd)
    append_log(model.subdir,"Finished: Calculating_phi_values")

def calculate_ddG_for_temperature(mutants,T,bounds,state_labels):
    """ Calculate ddG for"""
    if not os.path.exists("phi"):
        os.mkdir("phi")
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
    Q = np.loadtxt("Q.dat")
    state_indicator = np.zeros(len(Q),int)
    ## Assign every frame a state label. State indicator is integer 1-N for N states.
    for state_num in range(len(bounds)-1):
        instate = (Q > bounds[state_num]).astype(int)*(Q <= bounds[state_num+1]).astype(int)
        state_indicator[instate == 1] = state_num+1
    if any(state_indicator == 0):
        print "ERROR! There are unassigned frames!"
        print sum((state_indicator == 0).astype(int)), " unassigned frames out of ", len(Q)
        print " Exiting"
        raise SystemExit

    ## Boolean arrays that indicate which state each frame is in.
    ## States are defined by their boundaries along coordinate Q.
    U  = ((Q > bounds[1]).astype(int)*(Q < bounds[2]).astype(int)).astype(bool)
    TS = ((Q > bounds[3]).astype(int)*(Q < bounds[4]).astype(int)).astype(bool)
    N  = ((Q > bounds[5]).astype(int)*(Q < bounds[6]).astype(int)).astype(bool)
    Nframes  = float(sum(N.astype(int)))
    Uframes  = float(sum(U.astype(int)))
    TSframes = float(sum(TS.astype(int)))
    
    ## Compute deltaG for each state. Then DeltaDelta G with respect to the
    ## first state (assumed to be the unfolded/denatured state).
    ## Units of kT.
    dG = [[],[],[]]
    ddG = [[],[]]
    phi = []
    for k in range(len(mutants)):
        mut = mutants[k]
        print "    calculating ddG for ",mut
        ## Load dH_k for mutation. Eq. (2) minus Eq. (1) from reference (1).
        dH = np.loadtxt("dH_"+mut+".dat")

        ## Free energy perturbation formula. Equation (4) in reference (1).
        dG_U  = -np.log(np.sum(np.exp(-beta*dH[U]))/Uframes)
        dG_TS = -np.log(np.sum(np.exp(-beta*dH[TS]))/TSframes)
        dG_N  = -np.log(np.sum(np.exp(-beta*dH[N]))/Nframes)

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
    
    save_phi_values(mutants,"Q",state_labels,dG,ddG,phi)

def get_mutant_dH(path,mutants):
    """ Load the mutant energy perturbations dH_<mut>.dat

    Soon to be DEPRECATED 7-6-2014
    """

    i = 0
    for mut in mutants:
        temp = np.loadtxt(path+"/dH_"+mut+".dat")
        print "    Loading:",mut
        if i == 0:
            dH = np.zeros((len(mutants),len(temp)),float)
            dH[i,:] = temp
        else:
            dH[i,:] = temp
        i += 1
    
    return dH

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
    """ Calculates the normalized interaction betwen nonbonded pairs."""
    print "  Calculating Qij..."
    qij = model.nonbond_interaction(r,sig,delta)
    return qij

def get_state_bounds():
    """ Bounds for each state. Bounds are bin edges along Q. """
    if os.path.exists("whamQ/state_bounds.txt"):
        statefile = open("whamQ/state_bounds.txt","r").readlines()
    elif os.path.exists("state_bounds.txt"):
        statefile = open("state_bounds.txt","r").readlines()
    else:
        print "ERROR!"
        print "  Please create state_bounds.txt or whamQ/state_bounds.txt"
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

def get_Tf_choice():
    if os.path.exists("whamQ/Tf.txt"):
        Tf = open("whamQ/Tf.txt").read().rstrip("\n")
        print "  Using temp ",Tf
    elif os.path.exists("Tf.txt"):
        Tf = open("Tf.txt").read().rstrip("\n")
        print "  Using temp ",Tf
    else:
        print "ERROR!"
        print "  Please create Tf.txt or whamQ/Tf.txt with your choice to do mutations at."
        print "  Exiting"
        raise SystemExit
    if len(Tf.split("_")) > 1:
        Tnum = int(Tf.split("_")[1])
        Tf = float(Tf.split("_")[0])
    else:
        Tf = float(Tf)
        Tnum = 1
    return Tf,Tnum

def load_beadbead(subdir):
    """ Should be DEPRECATED """
    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt(subdir+"/BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    keep_interactions = np.zeros(len(interaction_numbers),int)
    for i in range(len(interaction_numbers)):
        if interaction_numbers[i] in ["ds","ss"]:
            pass
        else:
            keep_interactions[i] = int(interaction_numbers[i])

    return beadbead,keep_interactions
        
def load_eps_delta_sig_traj(subdir):
    """ Load in the info from the BeadBead.dat file. Sig_ij, eps_ij, delta_ij and
        index pairs. This information is constant for a trajectory. Filter all fields
        to keep only interactions with nonzero interaction type.

        In calculating the mutations only modify parameters that have interaction_type
        in the BeadBead.dat =/= [0,ds,ss]. 

        SOON TO BE DEPRECATED

    """
    print "  Loading BeadBead.dat"
    beadbead = np.loadtxt(subdir+"/BeadBead.dat",dtype=str) 
    sigij = beadbead[:,5].astype(float)
    epsij = beadbead[:,6].astype(float)
    deltaij = beadbead[:,7].astype(float)
    interaction_numbers = beadbead[:,4].astype(str)
    pairs = beadbead[:,:2].astype(int) 
    pairs -= np.ones(pairs.shape,int)

    keep_interactions = np.zeros(len(interaction_numbers),int)
    for i in range(len(interaction_numbers)):
        if interaction_numbers[i] in ["ds","ss"]:
            pass
        else:
            keep_interactions[i] = int(interaction_numbers[i])

    #print keep_interactions != 0       ## DEBUGGING
    #print sum((keep_interactions != 0).astype(int))      ## DEBUGGING
    sigij = sigij[keep_interactions != 0]
    epsij = epsij[keep_interactions != 0]
    deltaij = deltaij[keep_interactions != 0]
    pairs = pairs[keep_interactions != 0]

    print "  Only modifying ",sum((keep_interactions != 0).astype(int)), " parameters out of ", len(keep_interactions)
    ## Use mdtraj to compute the distances between pairs.
    print "  Loading traj.xtc with mdtraj..."
    traj = md.load(subdir+"/traj.xtc",top=subdir+"/Native.pdb")
    print "  Computing distances with mdtraj..."
    traj_dist = md.compute_distances(traj,pairs)

    return sigij,epsij,deltaij,interaction_numbers,keep_interactions,pairs,traj,traj_dist

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
    """ To Do: Make into a command line utility.
    """

    #parser = argparse.ArgumentParser(description='Calculate .')
    #parser.add_argument('--subdir', type=str, required=True, help='Directory.')
    #parser.add_argument('--calc_dH', action='store_true', help='Calculate dH for mutants.')
    #args = parser.parse_args()
    #pdb = args.subdir+".pdb"
    
    def dummy(this,that):
        pass

    pdb = "r15.pdb"
    model = models.SmogCalpha.SmogCalpha(pdb)
    Fij,Fij_pairs,Fij_conts,eps,deltas,sigmas,rij = calculate_dH_for_mutants(model,dummy)
    #calculate_phi_values(model,dummy)
