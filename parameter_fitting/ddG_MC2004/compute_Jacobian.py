''' Submodule to calculate the simulation phi values

Description:

    Submodule to calculate the energetic perturbation from each
mutation and the results DeltaDetla G's for simulation.


References:

(1) Matysiak, S.; Clementi, C. Optimal Combination of Theory and Experiment for
the Characterization of the Protein Folding Landscape of S6: How Far Can a
Minimalist model Go? J. Mol. Biol. 2004, 343, 235-248.
'''

import os
import argparse
import time
import numpy as np

import mdtraj as md

from mutatepdbs import get_core_mutations, get_scanning_mutations, get_exp_ddG

#import project_tools.parameter_fitting.util.util as util
from project_tools.parameter_fitting.util.util import *


global GAS_CONSTANT_KJ_MOL
GAS_CONSTANT_KJ_MOL = 0.0083144621

def get_dHk(model,rij,Fij_conts,Fij):
    ''' Get perturbed potential energy '''
    dHk = np.zeros(rij.shape[0],float)
    for i in range(len(Fij_conts)):       ## loop over number of parameters
        cont_idx = Fij_conts[i]
        dHk = dHk - Fij[i]*model.pairwise_strengths[cont_idx]*model.pairwise_potentials[cont_idx](rij[:,cont_idx])
    return dHk

def get_Vp_plus_Vpk(model,Vp,rij,Fij_conts,Fij):
    ''' Get perturbed potential energy '''
    Vp_plus_Vpk = np.array(Vp,copy=True)
    for i in range(len(Fij_conts)):       ## loop over number of parameters
        cont_idx = Fij_conts[i]
        param_idx = model.pairwise_param_assignment[cont_idx]
        Vp_plus_Vpk[:,param_idx] = Vp_plus_Vpk[:,param_idx] - \
                    Fij[i]*model.pairwise_potentials[cont_idx](rij[:,cont_idx])
    return Vp_plus_Vpk

def get_target_feature(model):
    ''' Get target features '''
    name = model.subdir
    iteration = model.iteration

    cwd = os.getcwd()
    os.chdir("%s/mutants" % name)
    target_feature, target_feature_err = get_exp_ddG()
    os.chdir(cwd)

    return target_feature, target_feature_err

def calculate_average_Jacobian(model,scanning_only=False,scanfij=0.5,saveas="Q_phi.dat",test=False):
    ''' Calculate the average feature vector (ddG's) and Jacobian '''
    
    name = model.subdir
    iteration = model.iteration

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    os.chdir("%s/mutants" % name)
    ## Get list of mutations and fraction of native contacts deleted for 
    ## each mutation.
    mutants_core = get_core_mutations()
    Fij_core, Fij_pairs_core, Fij_conts_core = get_mutant_fij(model,mutants_core)
    mutants_scanning = get_scanning_mutations()
    Fij_scanning, Fij_pairs_scanning, Fij_conts_scanning = get_mutant_fij_scanning(model,mutants_scanning,fij=scanfij)

    if scanning_only:
        mutants = mutants_scanning
        Fij = Fij_scanning
        Fij_pairs = Fij_pairs_scanning
        Fij_conts = Fij_conts_scanning
        saveas = "scan_%.2f.dat" % scanfij
    else:
        mutants = mutants_core + mutants_scanning
        Fij = Fij_core + Fij_scanning
        Fij_pairs = Fij_pairs_core + Fij_pairs_scanning
        Fij_conts = Fij_conts_core + Fij_conts_scanning

    os.chdir(sub)
    temperatures = [ x.split('_')[0] for x in open("long_temps_last","r").readlines() ] 
    directories = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ] 

    bounds, state_labels = get_state_bounds()
    bounds = [0] + bounds + [model.n_contacts]

    ## Loop over temperatures in iteration subdir. Calculate ddG vector and 
    ## Jacobian for each directory indpendently then save. Save the average
    ## feature vector and Jacobian in the iteration/newton directory.
    sim_feature_all = []
    Jacobian_all = []
    lasttime = time.time()
    for n in range(len(directories)):
        T = temperatures[n]
        dir = directories[n]
        beta = 1./(GAS_CONSTANT_KJ_MOL*float(T))
        print "  Calculating Jacobian for iteration_%d/%s" % (model.iteration,dir)
        os.chdir(dir)
        sim_feature, Jacobian = compute_Jacobian_for_directory(model,beta,mutants,Fij,Fij_pairs,Fij_conts,bounds,state_labels,saveas=saveas,test=test)
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

def compute_Jacobian_for_directory(model,beta,mutants,Fij,Fij_pairs,Fij_conts,bounds,state_labels,saveas="Q_phi.dat",test=False):
    ''' Calculates the feature vector (ddG's) and Jacobian for one directory '''
    ## Get trajectory, state indicators, contact energy
    traj,rij,Vp = get_rij_Vp(model)

    Q = np.loadtxt("Q.dat")
    #U,TS,N,Uframes,TSframes,Nframes = util.get_state_indicators(Q,bounds)
    U,TS,N,Uframes,TSframes,Nframes = get_state_indicators(Q,bounds)

    ## Average dimensionless potential energy for each state.
    Vp_U  = sum(Vp[U,:])/Uframes
    Vp_TS = sum(Vp[TS,:])/TSframes
    Vp_N  = sum(Vp[N,:])/Nframes

    ## Compute deltaG for each state. Then DeltaDelta G with respect to the
    ## first state (assumed to be the unfolded/denatured state).
    ## Units of kT.
    dG = [[],[],[]]
    ddG = [[],[]]
    phi = []

    ## Initialize Jacobian
    #Jacobian = np.zeros((2*len(mutants),model.n_contacts),float)
    Jacobian = np.zeros((2*len(mutants),model.n_model_param),float)
    sim_feature = np.zeros(2*len(mutants),float)

    avg_rowtime = []
    lasttime = time.time()
    for k in range(len(mutants)):
        mut = mutants[k]
        ## Compute energy perturbation
        dHk = get_dHk(model,rij,Fij_conts[k],Fij[k])
        np.savetxt("dH_%s.dat" % mut,dHk)

        ## Get perturbed dimensionless potential energy.
        Vp_plus_Vpk = get_Vp_plus_Vpk(model,Vp,rij,Fij_conts[k],Fij[k])

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

        ## simulation feature 
        sim_feature[k] = ddG_dagg
        sim_feature[k + len(mutants)] = ddG_stab

        ## Thermal averages for matrix equation (9).
        expdHk_U  = sum(np.exp(-beta*dHk[U]))/Uframes
        expdHk_TS = sum(np.exp(-beta*dHk[TS]))/TSframes
        expdHk_N  = sum(np.exp(-beta*dHk[N]))/Nframes

        ## Loop over number of model parameters. Could be done vectorized but idk if thats faster
        for p in range(model.n_model_param):
            Vp_Vpk_expdHk_U  = sum(Vp_plus_Vpk[U,p]*np.exp(-beta*dHk[U]))/Uframes
            Vp_Vpk_expdHk_TS = sum(Vp_plus_Vpk[TS,p]*np.exp(-beta*dHk[TS]))/TSframes
            Vp_Vpk_expdHk_N  = sum(Vp_plus_Vpk[N,p]*np.exp(-beta*dHk[N]))/Nframes

            Jacobian[k,p] = beta*(((Vp_Vpk_expdHk_TS/expdHk_TS) - Vp_TS[p]) - ((Vp_Vpk_expdHk_U/expdHk_U) - Vp_U[p]))

            Jacobian[k + len(mutants),p] = beta*(((Vp_Vpk_expdHk_N/expdHk_N) - Vp_N[p]) - ((Vp_Vpk_expdHk_U/expdHk_U) - Vp_U[p]))

        thistime = time.time()
        dt = thistime - lasttime
        avg_rowtime.append(dt)
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

def get_mutant_fij(model,mutants):
    ''' Load in the fraction of contact loss for each mutation.

    Description:

        Since the fraction of contacts lost matrix f^k_ij is sparse only load
    in the of nonzero elements and their indices. Determine the contact indices
    for nonzero entries. Let only allow fij's for contacts. 
    '''
    Fij_pairs = []
    Fij_conts = []
    Fij = []
    for mut in mutants:
        if model.nonnative:
            fij_temp = np.loadtxt("fij_%s.dat" % mut)
        else:
            fij_temp = model.Qref*np.loadtxt("fij_%s.dat" % mut)
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

def get_mutant_fij_scanning(model, mutants, fij=0.5):
    '''Estimate the local contact fraction loss for scanning mutations.
    
    The ddGs for scanning mutations could be affected by a multiplicity of factors that exceed the mere contacts lost (e.g. effect 
    of solvent, interaction with helix dipole, interaction with charged residues, loss of hydrogen bonds, among others). As a result,
    the exact determination varies according to which factor(s) weigh the most. Given that this is a coarse-grained model, as a 
    first approach we will estimate an average value of fij = 0.5. This is as said before arbitrary, but should give us an estimate
    of how the local epsilons vary according to the value of the ddGs involved.

    Also, on a first approach only the [i, i+4] contacts will be affected by the said value of fij.

    '''
    
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

            for j in range(model.n_contacts):
                if (model.contacts[j,0] == mut_res_number) and (model.contacts[j,1] == (mut_res_number+4)):
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
    ''' Save the calculated dG, ddG, and phi values for states'''

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
    return header_string+"\n"+data_string

if __name__ == "__main__":
    import model_builder as mdb

    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--name', type=str, required=True, help='name.')
    parser.add_argument('--iteration', type=int, required=True, help='iteration.')
    args = parser.parse_args()
    
    name = args.name
    iteration= args.iteration

    contacts = np.loadtxt("%s/contacts.dat" % name)
    pdb = "%s.pdb" % name
    defaults = True
    model = mdb.models.SmogCalpha.SmogCalpha(pdb=pdb,contacts=contacts,defaults=defaults,iteration=iteration)
    sim_feature_avg, sim_feature_err, Jacobian_avg, Jacobian_err = calculate_average_Jacobian(model,test=True)
