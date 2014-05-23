"""
Feb 14 2014
Alexander Kluber

    This submodule uses information in a mutations.dat file and the wild type
structure wt.pdb to create mutated pdbs and resulting contact maps. The contact
maps are then used to calculate the fraction of heavy-atom contact loss between
residues i and j for mutation k.
    Currently works! 2-14-14

MODELLER broke my Numpy build :( ===> FIXED by adding library path to LD_LIBRARY_PATH.
Follow instructions at:
https://docs.rice.edu/confluence/display/ITDIY/How+to+use+BLAS+and+LAPACK+libraries

"""

import numpy as np
import subprocess as sb
import os
import shutil

import mutatepdbs

def prep_mutants(System,append_log):
    """ Creates a mutated pdb for every mutant."""

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
    mutatepdbs.make_all_mutations() 
    print "  Calculating fraction of contact loss fij..."
    calculate_contacts_lost_for_mutants()

    os.chdir(cwd)

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

def get_res_res_conts(name):
    """ Get number of residue-residue heavy atom contacts from 
        all-atom contact map output from Shadowmap."""
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
    C90 = get_res_res_conts(name+".cutoff")
    diff = (Cwt - C90)
    print "    Number of contacts lost for ",name,sum(sum(diff))
    Cwt[ Cwt < 1 ] = 1.
    diff /= Cwt
    np.savetxt("fij_"+name+".dat",diff,fmt="%.5f",delimiter=" ")

    #plt.subplot(1,1,1,aspect=1)
    #plt.pcolor(diff,cmap=plt.cm.binary)
    #plt.xlabel("Residue $i$")
    #plt.ylabel("Residue $j$")
    #cbar = plt.colorbar()
    #cbar.set_label("Fraction contact loss $f_{ij}$")
    #plt.title("Fraction contact loss F90A")
    #plt.savefig("fij_wt_F90A.pdf")
    #plt.show()

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

    ## TO DO: Implement new mutations data file format.

    modelname = 'wt.pdb'
    mutation_data = open("mutations.dat","r").readlines()[1:]
    mut_indx = [ mutation_data[i].split()[0] for i in range(len(mutation_data)) ]
    wt_res = [ mutation_data[i].split()[1] for i in range(len(mutation_data)) ]
    mut_res = [ mutation_data[i].split()[2] for i in range(len(mutation_data)) ]
    
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

if __name__ == '__main__':
    calculate_contacts_lost_for_mutants()
