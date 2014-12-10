import numpy as np
import glob 


all_fij = glob.glob("fij*dat")

for i in range(len(all_fij)):
    print " writing mutation %d" % i
    fij_temp = np.loadtxt(all_fij[i])
    indices = np.nonzero(fij_temp)
    new_fij = np.ones(fij_temp.shape,float)
    new_fij[indices] = fij_temp[indices]
    np.savetxt(all_fij[i],new_fij)
