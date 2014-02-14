import numpy as np
import os
import subprocess as sb
import numpy as np


'''
Feb 5 2014
Alexander Kluber

    Starting to implement the mutations to the heterogeneous Go model.

MODELLER broke my Numpy build :( 

Follow instructions at:

https://docs.rice.edu/confluence/display/ITDIY/How+to+use+BLAS+and+LAPACK+libraries

'''



def get_shadow_pdb_atoms(name):
    ''' Parse the pdb file output by Shadow Jar.'''
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
    ''' Parses the output of shadow.jar to determine heavy atom contact
        map. Works. Requires the following files be present: name.contacts
        name.pdb
    '''
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
    ''' Get number of residue-residue heavy atom contacts from 
        all-atom contact map output from Shadowmap.'''
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
    ''' Calculate f^k_ij matrices for mutant.'''

    Cwt = get_res_res_conts("wt.cutoff")
    C90 = get_res_res_conts(name+".cutoff")
    diff = (Cwt - C90)
    print "Number of contacts lost for ",name,sum(sum(diff))
    Cwt[ Cwt < 1 ] = 1.
    diff /= Cwt
    #np.savetxt("fij_"+name+".dat",diff,fmt="%.5f",delimiter=" ")

    #plt.subplot(1,1,1,aspect=1)
    #plt.pcolor(diff,cmap=plt.cm.binary)
    #plt.xlabel("Residue $i$")
    #plt.ylabel("Residue $j$")
    #cbar = plt.colorbar()
    #cbar.set_label("Fraction contact loss $f_{ij}$")
    #plt.title("Fraction contact loss F90A")
    #plt.savefig("fij_wt_F90A.pdf")
    #plt.show()



if __name__ == '__main__':

    ## Start by reading the mutations file. Should be an attribute of System
    ## object.
    modelname = 'wt.pdb'
    mutation_data = open("mutations.txt","r").readlines()[1:]
    mut_indx = [ mutation_data[i].split()[0] for i in range(len(mutation_data)) ]
    wt_res = [ mutation_data[i].split()[1] for i in range(len(mutation_data)) ]
    mut_res = [ mutation_data[i].split()[2] for i in range(len(mutation_data)) ]
    
    cmd0 = 'cp /projects/cecilia/SCM.1.31.jar .'
    sb.call(cmd0,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
    cmd1 = 'echo -e "9\\n3\\n" | pdb2gmx -f wt.pdb -o wt.gro -p wt.top' 
    sb.call(cmd1,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
    cmd2 = 'java -jar SCM.1.31.jar -g wt.gro -t wt.top -o wt.cutoff.contacts -m cutoff -p wt.cutoff.wH.pdb'
    sb.call(cmd2,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))

    ## Use shadow map to create all-atom contact map. For each mutated pdb
    ## determine the 
    for i in range(len(mut_indx)):
        name = wt_res[i]+mut_indx[i]+mut_res[i]

        cmd2 = 'java -jar SCM.1.31.jar -g %s.gro -t %s.top -o %s.cutoff.contacts -m cutoff -p %s.cutoff.wH.pdb' % (name,name,name,name)
        #sb.call(cmd2,shell=True,stdout=open("cutoff.out","w"),stderr=open("cutoff.err","w"))
        
        print cmd1
        #calculate_fraction_contact_loss(name)
