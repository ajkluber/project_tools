import cython
import numpy as np
cimport numpy as np

@cython.boundscheck(False)
def count_transitions(int n_bins, int n_steps, np.ndarray[np.double_t, ndim=1] bins, np.ndarray[np.double_t, ndim=1] Q):

    cdef np.ndarray[np.double_t,
                    ndim=2,
                    negative_indices=False,
                    mode='c'] N = np.zeros((n_bins,n_bins))

    cdef int i,j,t

    for t in range(n_steps):
        for i in range(n_bins+1):
            if (bins[i] < Q[t]) and (bins[i+1] >= Q[t]):
                if t == 0:
                    j = i
                else:
                    N[i,j] += 1.
                    j = i
                break
    return N

@cython.boundscheck(False)
def calculate_M(int n_bins, np.ndarray[np.double_t, ndim=1] F, np.ndarray[np.double_t, ndim=1] D, np.double_t deltaQ):
    ''' Calculate matrix M from Bicout, Szabo 1998. 109'''

    cdef np.ndarray[np.double_t,
                    ndim=2,
                    negative_indices=False,
                    mode='c'] M = np.zeros((n_bins,n_bins))
    cdef int i

    for i in range(n_bins):
        if i == 0:
            M[i,i] = -((D[i+1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i+1])/2.)
            M[i,i+1] = np.sqrt(((D[i+1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i+1])/2.)*((D[i] + D[i+1])/(2.*deltaQ**2))*np.exp(-(F[i+1] - F[i])/2.))
            M[i+1,i] = M[i,i+1]
        elif i == (n_bins-1):
            M[i,i] = -((D[i-1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i-1])/2.)
            M[i,i-1] = np.sqrt(((D[i-1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i-1])/2.)*((D[i] + D[i-1])/(2.*deltaQ**2))*np.exp(-(F[i-1] - F[i])/2.))
            M[i-1,i] = M[i,i-1]
        else:
            M[i,i] =-((D[i+1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i+1])/2.) - ((D[i-1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i-1])/2.)
            M[i,i+1] = np.sqrt(((D[i+1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i+1])/2.)*((D[i] + D[i+1])/(2.*deltaQ**2))*np.exp(-(F[i+1] - F[i])/2.))
            M[i,i-1] = np.sqrt(((D[i-1] + D[i])/(2.*deltaQ**2))*np.exp(-(F[i] - F[i-1])/2.)*((D[i] + D[i-1])/(2.*deltaQ**2))*np.exp(-(F[i-1] - F[i])/2.))
            M[i+1,i] = M[i,i+1]
            M[i-1,i] = M[i,i-1]

    return M

@cython.boundscheck(False)
def calculate_propagator(int n_bins, np.ndarray[np.double_t, ndim=2] M, np.ndarray[np.double_t, ndim=1] F, np.double_t dt):

    cdef np.ndarray[np.double_t,
                    ndim=2,
                    negative_indices=False,
                    mode='c'] P = np.zeros((n_bins,n_bins))

    cdef np.ndarray[np.double_t,
                    ndim=2,
                    negative_indices=False,
                    mode='c'] vects = np.zeros((n_bins,n_bins))

    cdef np.ndarray[np.double_t,
                    ndim=1,
                    negative_indices=False,
                    mode='c'] vals = np.zeros(n_bins)

    cdef int i,j

    vals, vects = np.linalg.eig(M)

    for i in range(n_bins):
        vects[:,i] = vects[:,i]/np.linalg.norm(vects[:,i])

    for i in range(n_bins):
        for j in range(n_bins):
            for alpha in range(n_bins):
                P[i,j] += np.exp(-(F[j] - F[i])/2.)*vects[i,alpha]*vects[j,alpha]*np.exp(vals[alpha]*dt)

    return P


#@cython.boundscheck(False)
#def calculate_likelihood(int n_bins, np.ndarray[np.double_t, ndim=1] F, np.ndarray[np.double_t, ndim=1] D, np.double_t deltaQ ):
#
#    #cdef np.ndarray[np.double_t,
#    #                ndim=2,
#    #                negative_indices=False,
#    #                mode='c'] LogL = np.zeros((n_bins,n_bins))
#    cdef int i,j
#    cdef np.double_t LogL = 0.
#
#    ## 
#    for i in range(n_bins):
#        for j in range(n_bins):
#            LogL += N[i,j]*np P_m_n(F,vals,vects,i,j)
#
#    ## Apply a smoothening prior for D
#
#    print LogL
