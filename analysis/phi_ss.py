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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration


    phi = phi_prob.calculate_phi_values(name,iteration)

    element, bounds = get_sec_structure(name)

    phi_ss = []
    print "Element     Seq. Bounds  Phi"
    for i in range(len(element)):
        temp = np.mean(phi[bounds[i,0]:bounds[i,1]])
        phi_ss.append(temp)
        print "%-9s   %4d %4d   %.4f " % (element[i],bounds[i,0],bounds[i,1],temp)

