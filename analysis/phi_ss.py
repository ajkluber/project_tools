import os
import argparse
import numpy as np

import phi_prob


def get_sec_structure(name):
    
    element = []
    bounds = []
    for line in open("%s/secondary_structure.txt" % name,"r"):
        a = line.split()
        element.append(a[0])
        bounds.append([int(a[1]),int(a[2])])

    bounds = np.array(bounds)

    return element, bounds

def calculate_average_contact_phi_ss(name,iteration):

    phi = phi_prob.calculate_contact_phi_values(name,iteration)

    element, bounds = get_sec_structure(name)

    indxs = np.arange(len(phi))
    phi_ss = []
    print "Element     Seq. Bounds  Phi"
    output = "# Element     Seq. Bounds  Phi\n"
    for i in range(len(element)):
        ## only average over values that aren't exactly zero.
        residues = ((phi != 0.).astype(int)*(indxs >= (bounds[i,0] - 1)).astype(int)*(indxs <= (bounds[i,1] - 1)).astype(int)).astype(bool)
        #temp = np.mean(phi[bounds[i,0]:bounds[i,1]])
        temp = np.mean(residues)
        phi_ss.append(temp)
        print "%-9s   %4d %4d   %.4f " % (element[i],bounds[i,0],bounds[i,1],temp)
        output += "%-9s   %4d %4d   %.4f \n" % (element[i],bounds[i,0],bounds[i,1],temp)

    outfile = "%s/iteration_%d/contact_phi_ss" % (name,iteration)
    if not os.path.exists(outfile):
        open(outfile,"w").write(output)

    return phi_ss

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    calculate_average_phi_ss(name,iteration)

