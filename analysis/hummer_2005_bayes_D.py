''' Bayesian estimates of F(Q) and D(Q)



Ingredients:
  1. Reaction coordinate Q.
  2. Equilibrium trajectory Q(t)
  3. 

Procedure: 
    1. Assume uniform prior distribution of equilibrium probabili




'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

import szabo

def imshow_mask_zero(N):
    M = np.array(N,copy=True) 
    M[M == 0] = np.NaN
    jet = matplotlib.cm.jet
    temp = np.ma.array(M,mask=np.isnan(M))
    jet.set_bad('w',1.)
    plt.imshow(temp,cmap=jet,interpolation='none')
    plt.colorbar()
    plt.show()

n_bins = 30
timelag = 1 ## 500001 total frames

Qfull = np.loadtxt("Q.dat")
Q = Qfull[::timelag]
n_steps = len(Q)

bins = np.linspace(min(Q),max(Q)+1,num=n_bins+1)
deltaQ = bins[1] - bins[0]

print "Estimating diffusion model using n_steps: %d  n_bins: %d " % (n_steps,n_bins)
print "Count observed transitions between bins"
N = szabo.count_transitions(n_bins,n_steps,bins,Q)
#imshow_mask_zero(N)

print "Initializing guess for F(Q), D(Q)"
n,bins = np.histogram(Q,bins=bins)
F = -np.log(n) - min(-np.log(n))
D = np.ones(n_bins,float)

print "Calculating matrix M"
M = szabo.calculate_M(n_bins,F,D,deltaQ)
#imshow_mask_zero(M)

## Solve for the propagators using rate matrix or Szabo matrix diagonlization.
print "Calculating Propagator"
P = szabo.calculate_propagator(n_bins,M,F,float(timelag))

raise SystemExit
#print "Calculating Likelihood"
#vals, vects = np.linalg.eig(M)
#LogL = 0.
#for i in range(n_bins):
#    for j in range(n_bins):
#        LogL += N[i,j]*np.log(P[i,j])
#print LogL

## Monte Carlo optimize F(Q),D(Q) to maximize the log-likelihood function by
## exploring parameter space via Metropolis sampling procedure.
#n_iterations = 100
#for n in range(n_iterations):

    ## Generate candidate monte carlo move

    ## Calculate -lnL

    ## Accept or reject with Metropolis criterion















