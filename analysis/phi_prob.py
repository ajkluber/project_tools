import os
import argparse
import numpy as np

def calculate_contact_phi_values(name,iteration):
    ''' Calculate the percentage of native structure for each residue.'''

    pdb = open("%s/clean.pdb" % name,"r").readlines()

    n_residues = len(open("%s/Native.pdb" % name,"r").readlines()) - 1

    cwd = os.getcwd()
    sub =  "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)

    if not os.path.exists("contact_phi_res"):
        Tuse = open("long_temps_last","r").readlines()[0].rstrip("\n")
        contacts = np.loadtxt("%s/contacts.ndx" % Tuse,skiprows=1,dtype=int)

        #C = np.zeros((n_residues,n_residues),float)
        U  = np.loadtxt("cont_prob_%s.dat" % "U")
        TS = np.loadtxt("cont_prob_%s.dat" % "TS")
        N  = np.loadtxt("cont_prob_%s.dat" % "N")

        phi_res = np.zeros(n_residues,float)
        phi_contact = (TS - U)/(N - U)
        for i in range(n_residues):
            ## Average over contact formation for all contacts that residue i participates in.
            contacts_for_residue = ((contacts[:,0] == (i+1)).astype(int) + (contacts[:,1] == (i+1)).astype(int)).astype(bool)
            if np.any(contacts_for_residue):
                phi_res[i] = np.mean(phi_contact[contacts_for_residue])

        np.savetxt("contact_phi_res", phi_res)
    else:
        phi_res = np.loadtxt("contact_phi_res")

    if not os.path.exists("%s_contact_phi.pdb" % name):
        newpdb = ""
        for line in pdb:
            if line.startswith("END"):
                newpdb += "END"
                break
            else:
                resnum = int(line[22:26]) - 1
                newpdb += "%s     %5f\n" % (line[:55],phi_res[resnum])

        open("%s_contact_phi.pdb" % name,"w").write(newpdb)
    os.chdir(cwd)

    return phi_res

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    phi = calculate_contact_phi_values(name,iteration)






