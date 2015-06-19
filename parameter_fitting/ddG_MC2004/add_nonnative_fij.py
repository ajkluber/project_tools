import numpy as np
import glob 


all_fij = glob.glob("fij*dat")
nn_weight = 1.0
natives = np.loadtxt('../contact_map.dat')

#Does not make an exception for i:i+1, i+2 and i+3. Check that this is being taken care of

for i in range(len(all_fij)):
    print " writing mutation %d" % i
    fij_temp = np.loadtxt(all_fij[i])
    fij_natives = fij_temp + natives
    mutated_res = int(all_fij[i][5:-5]) - 1
    indices = np.nonzero(fij_natives)
    new_fij = np.zeros(fij_temp.shape)
    new_fij[mutated_res,:] = nn_weight * np.ones(fij_temp.shape[0])
    new_fij[:,mutated_res] = nn_weight * np.ones(fij_temp.shape[1])
    new_fij[indices] = fij_temp[indices]
    np.savetxt(all_fij[i],new_fij)
