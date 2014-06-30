""" Submodule to calculate the simulation phi values

Description:

    Submodule to calculate the energetic perturbation from each mutation and
the results DeltaDetla G's for simulation.


References:
"""

import numpy as np
import os
import argparse

import mdtraj as md

import model_builder.models as models
import model_builder.systems as systems

from mutatepdbs import get_core_mutations

global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def calculate_dH_for_mutants(model,append_log):
    """ Calculate Hamiltonian perturbation for each mutation dH_k """
    
    append_log(model.subdir,"Starting: Calculating_dH")

    cwd = os.getcwd()
    sub = cwd+"/"+model.subdir+"/Mut_"+str(model.Mut_iteration)
    ## Get the fraction of native contacts deleted for each mutation.
    os.chdir(model.subdir+"/mutants")
    mut_indx, wt_res, mut_res = get_core_mutations()
    mutants = [ wt_res[i]+mut_indx[i]+mut_res[i]  for i in range(len(mut_indx)) ]
    Fij, Fij_pairs, Fij_conts = get_mutant_fij(model,mutants)

    os.chdir(sub)

    #temperatures = [ x.split('_')[0] for x in open("T_array_last.txt","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ] 

    ## Loop over all directories in Mut subdirectory. For each trajectory
    ## and each mutation crunch dH.
    for dir in directories:
        os.chdir(dir)
        print "  Loading traj for: ", "Mut_"+str(model.Mut_iteration)+"/"+dir
        traj = md.load("traj.xtc",top="Native.pdb")
        for k in range(len(mutants)):
            if not os.path.exists("dH_"+mut+".dat"):
                print "    calc dH for ", mutants[k]
                mut = mutants[k]
                            
                ## Use mdtraj to compute the distances between pairs.
                rij = md.compute_distances(traj,Fij_pairs[k])
                
                ## Epsilons, deltas, and sigmas for relevant pairs for this mutation.
                eps = model.contact_epsilons[Fij_conts[k]]
                deltas = model.contact_deltas[Fij_conts[k]]
                sigmas = model.contact_sigmas[Fij_conts[k]]

                ## Calculate dH_k using distances and parameters. Save.
                Vij = eps*(5.*((sigmas/rij)**12) - 6.*((sigmas/rij)**10))
                dH_k = sum(Vij.T)
                np.savetxt("dH_"+mut+".dat",dH_k)
            else:
                pass
        os.chdir("..")
    os.chdir(cwd)
    append_log(model.subdir,"Finished: Calculating_dH")

def calculate_phi_values(Model,System,append_log,coord):
    """ Calculate the phi values for a trajectory. Requires only state 
        definitions and dH (energetic perturbation for each mutation).
    """
    
    append_log(System.subdir,"Starting: Calculating_phi_values")
    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/"+System.mutation_active_directory
    T = get_Tf_choice(sub)
    savedir = sub+"/"+T+"_agg"
    beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
    os.chdir(System.subdir)
    if not os.path.exists(savedir+"/phi"):
        os.mkdir(savedir+"/phi")

    os.chdir("mutants")
    mut_indx, wt_res, mut_res = get_core_mutations()
    os.chdir("..")
    mutants = [ wt_res[i]+mut_indx[i]+mut_res[i]  for i in range(len(mut_indx)) ]
    print "  Getting state bounds for coordinate:",coord
    bounds, states = get_state_bounds(savedir,coord)
    num_states = len(states)
    print "  Loading dH for mutants"
    dH = get_mutant_dH(savedir,mutants)

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
    append_log(System.subdir,"Finished: Calculating_phi_values")

def get_mutant_dH(path,mutants):
    """ Load the mutant energy perturbations dH_<mut>.dat"""

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

def get_exp_ddG():
    """ Get experimental ddG data from mutants/mutations.dat"""

    ddG_exp_all = np.loadtxt("mutants/mutations.dat",skiprows=1,usecols=(3,4))
     
    ddG_exp_TS_D = ddG_exp_all[:,0]
    ddG_exp_N_D = ddG_exp_all[:,1]

    return ddG_exp_TS_D, ddG_exp_N_D

def get_mutant_fij(model,mutants):
    """ Load in the fraction of contact loss for each mutation.

    Description:

        Since the fraction of contacts lost matrix f^k_ij is sparse only load
    in the of nonzero elements and their indices.
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
        for i in range(len(indices[0])):
            temppairs.append([indices[0][i],indices[1][i]])
            contact_num = list((model.contacts[:,0]-1 == indices[0][i]) & (model.contacts[:,1]-1 == indices[1][i])).index(True)
            tempconts.append(contact_num)
    
        Fij_conts.append(np.array(tempconts))
        Fij_pairs.append(temppairs)
    return Fij, Fij_pairs, Fij_conts

def get_Qij(Model,r,sig,delta,interaction_nums):
    """ Calculates the normalized interaction betwen nonbonded pairs."""
    print "  Calculating Qij..."
    qij = Model.nonbond_interaction(r,sig,delta)
    return qij

def get_state_bounds(path,coord):
    """ Get bounds for each state for specified coordinate. Return a list of boolean
        arrays that specifies if each frame is in the given state or not."""
    #print path+"/"+coord+"_states.txt" ## DEBUGGING
    #print open(path+"/"+coord+"_states.txt","r").read() ## DEBUGGING

    statefile = open(path+"/"+coord+"_states.txt","r").readlines()[1:]
    bounds = []
    for line in statefile:
        bounds.append([line.split()[0],float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    if coord in ["Q","Qnh","Qh"]:
        data = np.loadtxt(path+"/"+coord+".dat")
        data /= max(data)
    elif coord == "Nh":
        data = np.loadtxt(path+"/Nh.dat")
    elif coord == "Rg":
        dummy,data = np.loadtxt(path+"/radius_cropped.xvg",unpack=True)
    elif coord == "rmsd":
        dummy,data = np.loadtxt(path+"/rmsd.xvg",unpack=True)
    else:
        print "ERROR!"
        print "  No option for coordinate: ", coord
        print "  Exiting"
        raise SystemExit

    states = []
    for i in range(len(bounds)):
        print "  State: ", bounds[i][0], " is defined as between: ",bounds[i][2], bounds[i][3]
        states.append((bounds[i][2] <= data)*(data <= bounds[i][3]))

    return bounds,states

def get_Tf_choice(sub):
    if not os.path.exists(sub+"/Tf_choice.txt"):
        print "ERROR!"
        print "  Please create ",sub+"/Tf_choice.txt with your choice to do mutations at."
        print "  Exiting"
        raise SystemExit
    else:
        Tf_choice = open(sub+"/Tf_choice.txt").read().split()[0]
        print "  Calculating dH for temp ",Tf_choice
    return Tf_choice


def load_beadbead(subdir):

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

def save_phi_values(savedir,mutants,coord,bounds,dG,ddG,phi):
    """ Save the calculated dG, ddG, and phi values for states"""

    header_string = "# mut" 
    for i in range(len(bounds)):
        header_string += "     dG_"+bounds[i][0]+"(kT)"
    header_string += "    "
    for i in range(1,len(bounds)):
        header_string += " ddG_"+bounds[i][0]+"-"+bounds[0][0]+"(kT)"
    for i in range(1,len(bounds)-1):
        header_string += "   Phi_"+bounds[i][0]+"/"+bounds[-1][0]

    data_string = ''
    for j in range(len(mutants)):
        line = "%6s"%mutants[j]
        for i in range(len(dG)):
            line += "  %10.5f " % dG[i][j]
        for i in range(len(ddG)):
            line += "  %10.5f " % ddG[i][j]
        for i in range(len(phi)):
            line += "  %10.5f  " % phi[i][j]
        data_string += line+"\n"
    print "ddG and Phi values:"
    print header_string
    print data_string

    outputfile = open(savedir+"/phi/"+coord+"_phi.dat","w")
    outputfile.write(header_string+"\n"+data_string)
    outputfile.close()

if __name__ == "__main__":
    """ To Do: Make into a command line utility.
    """

    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--subdir', type=str, required=True, help='Directory.')
    parser.add_argument('--calc_dH', action='store_true', help='Calculate dH for mutants.')
    args.parse.parse_args()
    
    def dummy(this,that):
        pass
    model = models.SmogCalpha.SmogCalpha("r15.pdb")
    Fij, Fij_pairs, Fij_conts = calculate_dH_for_mutants(model,dummy)
