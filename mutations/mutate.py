import sys
import os

from modeller import *


'''
Feb 12 2014 
Alexander Kluber
    A hack of a MODELLER script for mutating residues that I found at the
following url:

http://salilab.org/modeller/wiki/Mutate%20model

Works like a charm out of the box.

** When on DaVinci (or any Rice server) must use the following calling style to
execute this:

modpy.sh python mutate.py

Where modpy.sh is a script installed by MODELLER that must be in your PATH.
Also, several paths need to be in your PYTHONPATH and LD_LIBRARY_PATH. See
MODELLER install script for specific paths.

Requires that the wild-type clean.pdb be present in the mutants/ directory
along with the mutations.txt file. 

Resulting pdbs will be named: <wt_res><mut_indx><mut_res>.pdb
For example mutating PHE (F) at position 90 in the wild-type structure to ALA
(A) would result in a pdb: F90A.pdb


'''

def residue_three_letter_code(rescode):
    '''Converting from three letter code to one letter FASTA code.'''
    residue_code = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN',
                    'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
                    'E': 'GLU', 'G': 'GLY', 'H': 'HIS',
                    'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
                    'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                    'S': 'SER', 'T': 'THR', 'W': 'TRP',
                    'Y': 'TYR', 'V': 'VAL'}
    return residue_code[rescode]

def make_all_mutations():
    ''' Read in mutational data from file. Parse the file into lists of mutation
        indices mut_indx, wt residue identity wt_res, and desired mutation 
        mut_res.
        
        MODELLER uses 3-letter amino acid code (e.g. ALA instead of A) whereas 
        the mutational data files will probably use single-letter code.
    '''
    modelname = 'wt'
    mutation_data = open("mutations.txt","r").readlines()[1:]
    mut_indx = [ mutation_data[i].split()[0] for i in range(len(mutation_data)) ]
    wt_res = [ mutation_data[i].split()[1] for i in range(len(mutation_data)) ]
    mut_res = [ mutation_data[i].split()[2] for i in range(len(mutation_data)) ]

    for i in range(len(mut_indx)):
        print "Performing mutation: %s%s%s" % (wt_res[i],mut_indx[i],mut_res[i])

        respos = mut_indx[i]
        restyp = residue_three_letter_code(mut_res[i])
        saveas = wt_res[i]+mut_indx[i]+mut_res[i]+".pdb"
        modeller_mutate_pdb(modelname,respos,restyp,saveas)


def modeller_mutate_pdb(modelname,respos,restyp,saveas,chain='A'):
    ''' Function to use MODELLER to mutate a pdb at index repos to residue
        restyp, then save the new pdb as saveas.

        Taken almost entirely as-is from a sample script in the MODELLER
        documentation. I deleted the final part that does energy minimization
        so that the backbone coordinates wouldn't be affected.'''

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

if __name__ == '__main__':
    make_all_mutations()
