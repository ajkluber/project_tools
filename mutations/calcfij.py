""" Calculate fraction of heavy atom contacts lost for mutation

"""

import os

import numpy as np
import matplotlib.pyplot as plt

from project_tools.mutations.mutatepdbs import get_all_core_mutations

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

def count_heavy_atom_contacts(pdb):
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
                if np.linalg.norm(atm_coords[i] - atm_coords[j]) <= 5.0:
                    C[res_nums[i]-1,res_nums[j]-1] += 1.0

    ## DEBUGGING
    #indices = np.nonzero(C)
    #for m in range(len(indices[0])):
    #    print "%5d %5d %5d" % (indices[0][m], indices[1][m], C[indices[0][m],indices[1][m]])

    return C, residues, atm_coords

if __name__ == "__main__":

    Cwt, wtresidues, wtcoords = count_heavy_atom_contacts("wt.pdb")

    #plt.pcolor(Cwt)
    #plt.colorbar()
    #plt.show()
    #raise SystemExit       ## DEBUGGING

    mutants = get_all_core_mutations()

    for k in range(len(mutants)):
    #for k in [11]:     ## DEUBUGGING
        mut = mutants[k]
        mutindx = int(mut[1:-1])
        print "Mutation: ",mut
        temp1 = np.nonzero(Cwt[mutindx-1,:])
        temp2 = np.nonzero(Cwt[:,mutindx-1])
        #print temp1,temp2
        #print Cwt[mutindx-1,temp1],Cwt[temp2,mutindx-1]

        Cmut, residues, coords = count_heavy_atom_contacts(mut+".pdb")

        print "%-8s%-10s%-4s%-5s%-7s%-5s" % ("Resi","Resj","Cwt","Cmut","diff","fij")
        Dmut = Cwt - Cmut
        indices = np.nonzero(Dmut)
        for m in range(len(indices[0])):
            a = indices[0][m]
            b = indices[1][m]
            #print wtresidues[a+1],wtresidues[b+1]
            #print residues[a+1], residues[b+1]
            fij = Dmut[a,b]/Cwt[a,b]
            print "%3s %-3d %3s %-3d  %3d  %3d  %3d  %9.5f" % (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij)

            #if a == mutindx-1:     ## DEBUGGING
            #    print "%3s*%-3d %3s %-3d  %3d  %3d  %3d  %9.5f" % (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij)
            #elif b == mutindx-1:
            #    print "%3s %-3d %3s*%-3d  %3d  %3d  %3d  %9.5f" % (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij)
            #else:
            #    print "%3s %-3d %3s %-3d  %3d  %3d  %3d  %9.5f" % (residues[a], a+1, residues[b], b+1, Cwt[a,b], Cmut[a,b], Dmut[a,b],fij)
        print "\n"

        #raise SystemExit       ## DEBUGGING
