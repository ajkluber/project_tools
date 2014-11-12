import cython 
import numpy as np
cimport numpy as np

@cython.boundscheck(False)
def generate_histo(int n_samples, int n_Ebins, np.ndarray[np.double_t, ndim=1] Ebins, np.ndarray[np.double_t, ndim=1] Ecumul,
                   int n_Xbins, np.ndarray[np.double_t, ndim=1] Xbins, np.ndarray[np.double_t, ndim=1] Xcumul):

    cdef np.ndarray[np.double_t,
                    ndim=2,
                    negative_indices=False,
                    mode='c'] gen_histo = np.zeros((n_samples,2))
    cdef np.ndarray[np.double_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] randoms = np.random.rand(n_samples)
    cdef int i,n
    
    for i in range(n_samples):
        for n in range(1,n_Ebins):
            if (randoms[i] > Ecumul[n-1]) and (randoms[i] <= Ecumul[n]):
                gen_histo[i,0] = 0.5*(Ebins[n-1] + Ebins[n])
                break
        for n in range(1,n_Xbins):
            if (randoms[i] > Xcumul[n-1]) and (randoms[i] <= Xcumul[n]):
                gen_histo[i,1] = 0.5*(Xbins[n-1] + Xbins[n])
                break

    return gen_histo
