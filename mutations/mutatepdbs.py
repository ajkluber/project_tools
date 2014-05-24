""" Mutate the pdbs with MODLLER

Description:

    This submodule uses information in a mutations.dat file and the wild type
structure wt.pdb to create mutated pdbs and resulting contact maps. The contact
maps are then used to calculate the fraction of heavy-atom contact loss between
residues i and j for mutation k.

MODELLER broke my Numpy build :( ===> FIXED by adding library path to LD_LIBRARY_PATH.
Follow instructions at:
https://docs.rice.edu/confluence/display/ITDIY/How+to+use+BLAS+and+LAPACK+libraries

"""
import subprocess as sb
import numpy as np
import shutil
import sys
import os

from modeller import *

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

def get_heavy_atom_contact_map(name):
    """ Parses the output of shadow.jar to determine heavy atom contact
        map. Works. Requires the following files be present: name.contacts
        name.pdb
    """
    atoms, num_heavy, heavy_indices, atms_per_res = get_shadow_pdb_atoms(name+".wH")
    wt_conts = np.loadtxt(name+".contacts",usecols=(1,3),dtype=int)
    C = np.zeros((num_heavy,num_heavy))

    for pair in wt_conts:
        atm1 = atoms[pair[0]-1]
        atm2 = atoms[pair[1]-1]

        if (atm1 != "BAD") and (atm2 != "BAD"):
            hvy1 = heavy_indices.index(pair[0])
            hvy2 = heavy_indices.index(pair[1])
            C[hvy1,hvy2] = 1
    return C, atms_per_res

def get_shadow_pdb_atoms(name):
    """ Parse the pdb file output by Shadow Jar."""
    prev_resid = 1
    num_heavy = 0
    atoms = []
    heavy_indices = []
    atms_per_res = []
    temp_num_atoms = 0
    for line in open(name+".pdb","r"):
        if line[:3] in ["TER","END"]:
            break
        else:
            atm = line[12:16].split()[0]
            atoms.append(atm)
            if atm != "BAD":
                index = int(line[6:11])
                heavy_indices.append(index)
                num_heavy += 1
                resid = int(line[22:26])
                if resid != prev_resid:
                    atms_per_res.append(temp_num_atoms)
                    temp_num_atoms = 1
                    prev_resid = resid
                else:
                    temp_num_atoms += 1

    atms_per_res.append(temp_num_atoms)
    return atoms,num_heavy,heavy_indices,atms_per_res

def get_core_mutations():
    """ Extract mutational data. Only return info for useable mutations """
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
    
    return mut_indx,wt_res,mut_res

def get_core_mutation_ddG():
    """ Extract mutational data. Only return info for useable mutations """
    mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
    useable_and_core = []
    for i in range(mutation_data.shape[0]):
        if (mutation_data[i,0] == "core") and (mutation_data[i,8] == "True") and (mutation_data[i,11] == "0"):
            useable_and_core.append(True)
        else:
            useable_and_core.append(False)
    useable_and_core = np.array(useable_and_core)
    
    ddG_N_D = mutation_data[(useable_and_core == True),4] 
    ddG_N_D_err = mutation_data[(useable_and_core == True),5] 
    ddG_TS_D = mutation_data[(useable_and_core == True),6] 
    ddG_TS_D_err = mutation_data[(useable_and_core == True),7] 
    
    return ddG_N_D,ddG_N_D_err,ddG_TS_D,ddG_TS_D_err

#def get_mutational_data():
#    """ Extract mutational data. Only return info for useable mutations """
#    mutation_data = np.loadtxt("calculated_ddG.dat",dtype=str)
#    useable = np.array([ bool(x) for x in mutation_data[:,8] ])
#    mut_indx = mutation_data[(useable == True),1] 
#    wt_res = mutation_data[(useable == True),2] 
#    mut_res = mutation_data[(useable == True),3] 
#    
#    return mutation_data,useable,mut_indx,wt_res,mut_res

def get_res_res_conts(name):
    """ Get number of residue-residue heavy atom contacts from 
        all-atom contact map output from Shadowmap.
    """
    C, atms_per_res = get_heavy_atom_contact_map(name)
    N = len(atms_per_res)
    C_res = np.zeros((N,N),float)
    for i in range(4,N):
        lindx = sum(atms_per_res[:i])
        rindx = sum(atms_per_res[:i+1])
        for j in range(i+4,N):
            bindx = sum(atms_per_res[:j])
            tindx = sum(atms_per_res[:j+1])
            res_res_conts = C[lindx:rindx,bindx:tindx]
            C_res[i,j] = float(sum(sum(res_res_conts)))
    return C_res

def calculate_fraction_contact_loss(name):
    """ Calculate f^k_ij matrices for mutant."""

    Cwt = get_res_res_conts("wt.cutoff")
    Cmut = get_res_res_conts(name+".cutoff")
    diff = (Cwt - Cmut)
    print "    Number of contacts lost for ",name,sum(sum(diff))
    Cwt[ Cwt < 1 ] = 1.
    diff /= Cwt
    np.savetxt("fij_"+name+".dat",diff,fmt="%.5f",delimiter=" ")

def calculate_contacts_from_pdb(name):
    """ Calls shadow map to calculate"""
    if os.path.exists(name+".gro") == False:
        cmd1 = 'echo -e "9\\n3\\n" | pdb2gmx -f %s.pdb -o %s.gro -p %s.top' % (name,name,name)
        sb.call(cmd1,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
    cmd2 = 'java -jar SCM.1.31.jar -g %s.gro -t %s.top -o %s.cutoff.contacts -m cutoff -p %s.cutoff.wH.pdb' % (name,name,name,name)
    sb.call(cmd2,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))

def calculate_contacts_lost_for_mutants():
    """ Calculates the fraction of heavy-atom contact loss between residues i and j
        for mutant k:  f^k_ij . Must be in mutants directory which holds wild-type 
        pdb wt.pdb, a file hold mutations information mutations.dat."""

    mut_indx,wt_res,mut_res = get_core_mutations()

    if os.path.exists("SCM.1.31.jar") == False:
        cmd0 = 'cp /projects/cecilia/SCM.1.31.jar .'
        sb.call(cmd0,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
    if os.path.exists("wt.cutoff.contacts") == False:
        if os.path.exists("wt.gro") == False:
            cmd1 = 'echo -e "9\\n3\\n" | pdb2gmx -f wt.pdb -o wt.gro -p wt.top' 
            sb.call(cmd1,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
        cmd2 = 'java -jar SCM.1.31.jar -g wt.gro -t wt.top -o wt.cutoff.contacts -m cutoff -p wt.cutoff.wH.pdb'
        sb.call(cmd2,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))

    ## Use shadow map to create all-atom contact map. For each mutated pdb
    ## determine the fraction of heavy-atom contacts lost .
    for i in range(len(mut_indx)):
        name = wt_res[i]+mut_indx[i]+mut_res[i]
        if os.path.exists(name+".cutoff.contacts") == False:
            print "    Calculating contacts for", name
            calculate_contacts_from_pdb(name)
        if os.path.exists("fij_"+name+".dat") == False:
            print "    Calculating fij for ", name
            calculate_fraction_contact_loss(name)


def make_all_mutations():
    """ Read in mutational data from file. Parse the file into lists of mutation
        indices mut_indx, wt residue identity wt_res, and desired mutation 
        mut_res.
        
        MODELLER uses 3-letter amino acid code (e.g. ALA instead of A) whereas 
        the mutational data files will probably use single-letter code.
    """

    modelname = 'wt'
    mut_indx,wt_res,mut_res = get_core_mutations()

    for i in range(len(mut_indx)):

        saveas = wt_res[i]+mut_indx[i]+mut_res[i]+".pdb"
        if not os.path.exists(saveas):

            print "    Performing mutation: %s%s%s" % (wt_res[i],mut_indx[i],mut_res[i])

            respos = mut_indx[i]
            restyp = residue_three_letter_code(mut_res[i])
            saveas = wt_res[i]+mut_indx[i]+mut_res[i]+".pdb"
            modeller_mutate_pdb(modelname,respos,restyp,saveas)
        else:
            print "    Skipping mutation: %s%s%s" % (wt_res[i],mut_indx[i],mut_res[i])

def modeller_mutate_pdb(modelname,respos,restyp,saveas,chain='A'):
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

    log.verbose()

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

def prepare_mutants(System,append_log):
    """ Creates a mutated pdb for every mutant."""

    append_log("Starting: Preparing_Mutants")
    cwd = os.getcwd()
    sub = cwd+"/"+System.subdir+"/mutants"
    mut_filename = System.subdir+"_calculated_ddG.dat"
    if not os.path.exists(sub):
        print "  Creating direcotory",sub
        os.mkdir(sub)
    if not os.path.exists("wt.pdb"):
        print "  Copying wt"
        shutil.copy(System.subdir+"/clean_noH.pdb","wt.pdb")
    if (not os.path.exists(sub+"/calculated_ddG.dat")):
        print "  Didn't find "+sub+"/calculated_ddG.dat. Looking for "+System.path+"/"+mut_filename
        if os.path.exists(System.path+"/"+mut_filename):
            shutil.copy(System.path+"/"+mut_filename,sub+"/calculated_ddG.dat")
        else:
            print "ERROR!"
            print "  Didn't find "+mut_filename+" in "+System.path
            print "  Prepare the mutational data using prepare_ddG.py"
            print "  Exiting."
            raise SystemExit

    os.chdir(sub)

    print "  Mutating pdbs with MODELLER..."
    make_all_mutations() 
    print "  Calculating fraction of contact loss fij..."
    calculate_contacts_lost_for_mutants()
    os.chdir(cwd)
    append_log("Finished: Preparing_Mutants")

def command_line_prepare_mutants(System,append_log):
    """ Creates a mutated pdb for every mutant."""

    if not os.path.exists("calculated_ddG.dat"):
        print "ERROR!"
        print "  The file calculated_ddG.dat must exist!"
        print "  exiting..."
        raise SystemExit

    print "  Mutating pdbs with MODELLER..."
    make_all_mutations() 
    print "  Calculating fraction of contact loss fij..."
    calculate_contacts_lost_for_mutants()

if __name__ == '__main__':
    command_line_prepare_mutants()
