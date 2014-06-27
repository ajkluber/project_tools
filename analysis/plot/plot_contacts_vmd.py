import numpy as np
import matplotlib.pyplot as plt

beadbead = np.loadtxt("BeadBead.dat",dtype=str)
pairs = beadbead[:,:2].astype(int)
pairs -= np.ones(pairs.shape)
epsij = beadbead[:,6].astype(float)
native = beadbead[:,4].astype(int)

pairs = pairs[ native != 0 ]
epsij = epsij[ native != 0 ]

Qref = np.loadtxt("../../Qref_cryst.dat")
N = len(Qref)
C = np.zeros((N,N),float)

for k in range(pairs.shape[0]):
    C[pairs[k,0],pairs[k,1]] = epsij[k]

C_per_res = sum(C)
C_per_res /= float(max(C_per_res))

pdbfile = open("../../clean.pdb","r").readlines()

newpdb = ""
for i in range(len(pdbfile)):
    if pdbfile[i].startswith("END"):
        break
    else:
        oldline = pdbfile[i][:-1]
        #resnum = int(oldline[22:26]) 
        resnum = int(oldline[22:26]) - 1
        newline = oldline + ("%5.2f" % C_per_res[resnum]) + "\n"
        newpdb += newline
        #print len(oldline)
        #print oldline

newpdb += "END"
open("residues_by_epsilons.pdb","w").write(newpdb)
    
