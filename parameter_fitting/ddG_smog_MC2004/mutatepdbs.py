""" Mutate the pdbs with MODELLER

Description:

    This submodule uses information in a mutations.dat file and the wild type
structure wt.pdb to create mutated pdbs and resulting contact maps. The contact
maps are then used to calculate the fraction of heavy-atom contact loss between
residues i and j for mutation k.

MODELLER broke my Numpy build :( ===> FIXED by adding library path to LD_LIBRARY_PATH.
Follow instructions at:
https://docs.rice.edu/confluence/display/ITDIY/How+to+use+BLAS+and+LAPACK+libraries

"""
import argparse
import numpy as np
import os

def residue_three_letter_code(rescode):
    """Converting from three letter code to one letter FASTA code."""
    residue_code = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN',
                    'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
                    'E': 'GLU', 'G': 'GLY', 'H': 'HIS',
                    'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
                    'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                    'S': 'SER', 'T': 'THR', 'W': 'TRP',
                    'Y': 'TYR', 'V': 'VAL'}
    return residue_code[rescode]

def get_all_core_mutations():
    """ Extract mutational data. Only return info for useable mutations """
    if os.path.exists("core.ddG"):
        data = np.loadtxt("core.ddG",skiprows=1,dtype=str)
        mutants = [ "%s%s%s" % (data[i,1],data[i,0],data[i,2]) for i in range(data.shape[0]) ]
    else:
        mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
        core_muts = []
        for i in range(mutation_data.shape[0]):
            if (mutation_data[i,0] == "core"):
                core_muts.append(True)
            else:
                core_muts.append(False)
        core_muts = np.array(core_muts)

        mut_indx = mutation_data[(core_muts == True),1] 
        wt_res = mutation_data[(core_muts == True),2] 
        mut_res = mutation_data[(core_muts == True),3] 
        mutants = [ wt_res[i]+mut_indx[i]+mut_res[i]  for i in range(len(mut_indx)) ]
    
    return mutants

def get_scanning_mutations():
    """ Return alanine-glycine scanning mutants."""
    if os.path.exists("scanning.ddG"):
        data = np.loadtxt("scanning.ddG",skiprows=1,dtype=str)
        use = np.where(data[:,8] == "0")[0]
        mutants = [ "%s%s%s" % (data[i,1],data[i,0],data[i,2]) for i in use ]
    elif os.path.exists("calculated_ddG.dat"):
        mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
        scan_muts = []
        for i in range(mutation_data.shape[0]):
            if (mutation_data[i,0] == "surf") and (mutation_data[i,8] == "True") and (mutation_data[i,11] == "0"):
                scan_muts.append(True)
            else:
                scan_muts.append(False)

        scan_muts = np.array(scan_muts)

        mut_indx = mutation_data[(scan_muts == True),1] 
        wt_res = mutation_data[(scan_muts == True),2] 
        mut_res = mutation_data[(scan_muts == True),3] 
        mutants = [ wt_res[i]+mut_indx[i]+mut_res[i]  for i in range(len(mut_indx)) ]
    else:
        mutants = []
    
    return mutants

def get_core_mutations():
    """ Extract mutational data. Only return info for useable mutations """
    if os.path.exists("core.ddG"):
        data = np.loadtxt("core.ddG",skiprows=1,dtype=str)
        use = np.where(data[:,8] == "0")[0]
        mutants = [ "%s%s%s" % (data[i,1],data[i,0],data[i,2]) for i in use ]
    else:
        mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
        useable_and_core = []
        for i in range(mutation_data.shape[0]):
            if (mutation_data[i,0] == "core") and (mutation_data[i,8] == "True") and (mutation_data[i,11] == "0"):
                useable_and_core.append(True)
            else:
                useable_and_core.append(False)
        useable_and_core = np.array(useable_and_core)

        mut_indx = mutation_data[(useable_and_core == True),1] 
        wt_res = mutation_data[(useable_and_core == True),2] 
        mut_res = mutation_data[(useable_and_core == True),3] 
        mutants = [ wt_res[i]+mut_indx[i]+mut_res[i]  for i in range(len(mut_indx)) ]
    
    return mutants

def get_core_mutation_ddG():
    """ Extract mutational data. Only return info for useable mutations """
    if os.path.exists("core.ddG"):
        data = np.loadtxt("core.ddG",skiprows=1,dtype=str)
        use = np.where(data[:,8] == "0")[0]
        mutants = [ "%s%s%s" % (data[i,1],data[i,0],data[i,2]) for i in use ]
        ddG_N_D = data[use,3].astype(float)
        ddG_N_D_err = data[use,4].astype(float)  
        ddG_TS_D = data[use,5].astype(float)      
        ddG_TS_D_err = data[use,6].astype(float)  
    else:
        mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
        useable_and_core = []
        for i in range(mutation_data.shape[0]):
            if (mutation_data[i,0] == "core") and (mutation_data[i,8] == "True") and (mutation_data[i,11] == "0"):
                useable_and_core.append(True)
            else:
                useable_and_core.append(False)
        useable_and_core = np.array(useable_and_core)
        
        ddG_N_D = mutation_data[(useable_and_core == True),4].astype(float)
        ddG_N_D_err = mutation_data[(useable_and_core == True),5].astype(float)
        ddG_TS_D = mutation_data[(useable_and_core == True),6].astype(float) 
        ddG_TS_D_err = mutation_data[(useable_and_core == True),7].astype(float) 
    
    return ddG_N_D,ddG_N_D_err,ddG_TS_D,ddG_TS_D_err

def get_exp_ddG():
    """ Return both surface and core ddG_exp that are useable """
    print "  Getting experimental ddG"
    if os.path.exists("core.ddG"):
        data = np.loadtxt("core.ddG",skiprows=1,dtype=str)
        use = np.where(data[:,8] == "0")[0]
        mutants = [ "%s%s%s" % (data[i,1],data[i,0],data[i,2]) for i in use ]
        ddG_N_D = data[use,3].astype(float)
        ddG_N_D_err = data[use,4].astype(float)  
        ddG_TS_D = data[use,5].astype(float)      
        ddG_TS_D_err = data[use,6].astype(float)  
    else:
        mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
        useable = []
        for i in range(mutation_data.shape[0]):
            #print "%s %s%s%s" % (mutation_data[i,0],mutation_data[i,2],mutation_data[i,1],mutation_data[i,3])
            if (mutation_data[i,8] == "True") and (mutation_data[i,11] == "0"):
                useable.append(True)
            else:
                useable.append(False)
        useable = np.array(useable)
        
        ddG_N_D = mutation_data[(useable == True),4].astype(float)
        ddG_N_D_err = mutation_data[(useable == True),5].astype(float)
        ddG_TS_D = mutation_data[(useable == True),6].astype(float) 
        ddG_TS_D_err = mutation_data[(useable == True),7].astype(float) 
        
    ddGexp = np.concatenate((ddG_TS_D,ddG_N_D),axis=0)
    ddGexp_err = np.concatenate((ddG_TS_D_err,ddG_N_D_err),axis=0) 
    
    return ddGexp, ddGexp_err

def get_sim_ddG(mutants,coord):
    """ Get the saved ddG from simulation that should have been computed already."""

    ddGsim_TS_D = ddG = np.zeros(len(mutants),float)
    ddGsim_N_D = ddG = np.zeros(len(mutants),float)
    ddG_sim_all = np.loadtxt("phi/"+coord+"_phi.dat",skiprows=1,usecols=(0,4,5),dtype=str)
    sim_muts = list(ddG_sim_all[:,0])
    for k in range(len(mutants)):
        try:
            temp_indx = sim_muts.index(mutants[k])
        except:
            print "ERROR!"
            print "  The ddG_simulation was not found for mutation ", mutants[k]
            print "  Double check that ddG's are there for all mutations used."
            print "  Exiting"
            raise SystemExit
        ddGsim_TS_D[k] = float(ddG_sim_all[temp_indx,1])
        ddGsim_N_D[k] = float(ddG_sim_all[temp_indx,2])
        #ddGsim_TS_D_err[k] = float(ddG_sim_all[temp_indx,1])
        #ddGsim_N_D_err[k] = float(ddG_sim_all[temp_indx,2])
        
    return ddGsim_TS_D, ddGsim_N_D

def get_AApdb_coords(pdb):
    """ Parse atom names, indices, coordinates; residue names, indices from all-atom pdb """

    ## Only interested in ATOM lines.
    atmlines = [ line.rstrip('\n') for line in open(pdb,'r').readlines() if line.startswith("ATOM") ] 

    atm_nums   = []
    atm_names  = []
    atm_coords = []
    res_nums   = []
    res_names  = []

    for line in atmlines:
        atom_type = line[12:16].strip()
        atom_type = atom_type.rstrip()

        ## Only keep heavy atoms
        if atom_type.startswith("H"):
            continue
        else:
            atm_nums.append(int(line[6:11]))
            atm_names.append(atom_type)
            atm_coords.append(np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])]))
            res_nums.append(int(line[22:26]))
            res_names.append(line[17:20])

    atm_coords = np.array(atm_coords)
    return atm_nums,atm_names,atm_coords,res_nums,res_names

def count_heavy_atom_contacts(pdb,cutoff=6.0):
    """ Calculate # of residue-residue heavy atom contacts. """

    atm_nums,atm_names,atm_coords,res_nums,res_names = get_AApdb_coords(pdb)
    n_atoms = len(res_nums)
    n_residues = max(res_nums)

    residues = []
    temp = -1
    for m in range(n_atoms):
        if res_nums[m] != temp:
            residues.append(res_names[m])
            temp = res_nums[m]
        else:
            continue
    #print residues      ## DEBUGGING

    C = np.zeros((n_residues,n_residues),float)
    for i in range(n_atoms):
        for j in range(i,n_atoms):
            if res_nums[j] <= (res_nums[i]+3):
                pass
            else:
                #print np.linalg.norm(atm_coords[i] - atm_coords[j])    ## DEBUGGING
                if np.linalg.norm(atm_coords[i] - atm_coords[j]) <= cutoff:
                    C[res_nums[i]-1,res_nums[j]-1] += 1.0
    ## DEBUGGING
    #indices = np.nonzero(C)
    #for m in range(len(indices[0])):
    #    print "%5d %5d %5d" % (indices[0][m], indices[1][m], C[indices[0][m],indices[1][m]])

    return C, residues, atm_coords

def calculate_contacts_lost_for_mutants(avgflag):
    """ Calculates the fraction of heavy-atom contact loss between residues i and j
        for mutant k:  f^k_ij . Must be in mutants directory which holds wild-type 
        pdb wt.pdb."""

    Cwt, wtresidues, wtcoords = count_heavy_atom_contacts("wt.pdb")
    if not os.path.exists("Cwt.dat"):
        np.savetxt("Cwt.dat",Cwt)

    ## Grab native contact map
    Cmap = np.zeros(Cwt.shape,int)
    contacts = np.loadtxt("contacts",dtype=int) 
    if contacts.shape[1] == 4:
        contacts = contacts[:,0:3:2] ## Is that right?
    for pair in contacts:
        Cmap[pair[0]-1,pair[1]-1] = 1

    modelname = 'wt'
    mutants = get_all_core_mutations()
    log_string = "Core Mutations:\n"
    for k in range(len(mutants)):
        mut = mutants[k]
        print "    Calculating fij for: %s" % mut
        respos = mut[1:-1]
        restyp = residue_three_letter_code(mut[-1])
        saveas = "%s.pdb" % mut
        modeller_mutate_pdb(modelname,respos,restyp,saveas)

        mutindx = int(mut[1:-1])

        Cmut, residues, coords = count_heavy_atom_contacts(mut+".pdb")
        tempstring = "Mutation: %s\n" % mut
        tempstring += "%-8s%-10s%-4s%-5s%-7s%-5s\n" % ("Resi","Resj","Cwt","Cmut","diff","fij")
        Dmut = (Cwt - Cmut)*Cmap
        indices = np.nonzero(Dmut)
        Fij = np.zeros(Cwt.shape)
        nonzerofij = [ Dmut[indices[0][m],indices[1][m]]/Cwt[indices[0][m],indices[1][m]] for m in range(len(indices[0])) \
                            if Cwt[indices[0][m],indices[1][m]] != 0 ]
        avgfij =  np.mean(nonzerofij)
        for m in range(len(indices[0])):
            a = indices[0][m]
            b = indices[1][m]
            if Cwt[a,b] == 0.:
                continue
            else:
                if avgflag:
                    fij = Dmut[a,b]/Cwt[a,b]
                    Fij[a,b] = avgfij
                    tempstring += "%3s %-3d %3s %-3d  %3d  %3d  %3d  %9.5f %9.5f\n" % \
                        (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij,avgfij)
                else:
                    fij = Dmut[a,b]/Cwt[a,b]
                    Fij[a,b] = fij
                    tempstring += "%3s %-3d %3s %-3d  %3d  %3d  %3d  %9.5f\n" % \
                        (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij)
        np.savetxt("fij_%s.dat" % mut,Fij,fmt="%.6f",delimiter=" ")
        log_string += tempstring+"\n"

    log_string += "\nAlanine-Glycine Mutations:\n"
    scanmuts = get_scanning_mutations()
    for k in range(len(scanmuts)): 
        mut = scanmuts[k]
        mutindx = int(mut[1:-1])

    open("mutatepdbs.log","w").write(log_string)

def modeller_mutate_pdb(modelname,respos,restyp,saveas,chain=''):  #Replaced chain='A' by chain='' FY 11-NOV-2015
    """ Use MODELLER to mutate the pdb modelname at respos to restyp then save as saveas.

        Taken almost entirely as-is from a sample script in the MODELLER
    documentation. I deleted the final part that does energy minimization
    so that the backbone coordinates wouldn't be affected.

    A hack of a MODELLER script for mutating residues that I found at the
    following url:

    http://salilab.org/modeller/wiki/Mutate%20model

    Works like a charm out of the box.

    ** When on DaVinci (or any Rice server) must use the following calling style to
    execute this:

    modpy.sh python2.6 mutate.py

    Where modpy.sh is a script installed by MODELLER that must be in your PATH.
    Also, several paths need to be in your PYTHONPATH and LD_LIBRARY_PATH. See
    MODELLER install script for specific paths.

    Requires that the wild-type clean.pdb be present in the mutants/ directory
    along with the mutations.dat file. 

    Resulting pdbs will be named: <wt_res><mut_indx><mut_res>.pdb
    For example mutating PHE (F) at position 90 in the wild-type structure to ALA
    (A) would result in a pdb: F90A.pdb
    """

    log.none()

    # Set a different value for rand_seed to get a different final model
    env = environ(rand_seed=-49837)

    env.io.hetatm = True
    #soft sphere potential
    env.edat.dynamic_sphere=False
    #lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = model(env, file=modelname)
    ali = alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    #set up the mutate residue selection segment
    s = selection(mdl1.chains[chain].residues[respos])

    #perform the mutate residue operation
    s.mutate(residue_type=restyp)


    #get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])


    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    #here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    #yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = model(env, file=modelname)

    #required to do a transfer_res_numb
    #ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    #transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

    #It is usually necessary to write the mutated sequence out and read it in
    #before proceeding, because not all sequence related information about MODEL
    #is changed by this command (e.g., internal coordinates, charges, and atom
    #types and radii are not updated).

    mdl1.write(file=saveas)

if __name__ == '__main__':
    """ Creates a mutated pdb for every mutant."""

    
    parser = argparse.ArgumentParser(description='Calculate .')
    parser.add_argument('--avg',type=bool, default=False, help='Use avgfij.')
    args = parser.parse_args()
    
    avgflag = args.avg


    from modeller import *
    if (not os.path.exists("calculated_ddG.dat")) and (not os.path.exists("core.ddG")):
        raise IOError("calculated_ddG.dat or core.ddG must exist!")

    print "  Calculating fraction of contact loss fij."
    if avgflag:
        print "  Using average fij per mutation"
    calculate_contacts_lost_for_mutants(avgflag)
    print "  See mutatepdbs.log for summary of contact loss per mutant."
