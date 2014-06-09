""" Estimate convergence of standard error for simulation observables.


Description:

    There is no global criteria to prove that a simulations -has- converged,
but there are many ways to tell that a simulations -has not- converged. This
submodule plots some criteria discussed in reference (1).

    For example, the observables of a system in thermal equilibrium have 
stationary expectation values <A_i>. 

References:

(1) Grossfield, A.; Zuckerman, D. Quantifying uncertainty and sampling quality
in biomolecular simulations. Annu. Rep. Comput. Chem. 2009, 5, 23-48.
"""


import numpy as np
import matplotlib.pyplot as plt

def plot_BSE(Q):
    BSE = []
    for n in range(1,len(Q)-1):
        if (n % 5000) == 0:
            print n
        M = len(Q)/n
        if M <= 10:
            break
        temp = [ np.mean(Q[k*n:(k+1)*n]) for k in range(M) ]
        BSE.append(np.std(temp)/np.sqrt(float(M)))

    plt.plot(BSE)
    plt.show()


def plot_avg_versus_window_size():
    T = "148.76"
    #for i in range(1,6):
    for i in [1]:

        #Q = np.loadtxt(T+"_"+str(i)+"/Q.dat")
        Q = np.loadtxt(T+"_agg/Q.dat")
        #t,Q = np.loadtxt(T+"_"+str(i)+"/radius_cropped.xvg",unpack=True)
        #t,Q = np.loadtxt(T+"_"+str(i)+"/rmsd.xvg",unpack=True)
        #t,Q = np.loadtxt(T+"_"+str(i)+"/phis.xvg",unpack=True,usecols=(3,4))

        #plt.plot(Q)
        #plt.show()

        Avg = np.array([ np.mean(Q[:m]) for m in range(1,len(Q)) ])
        Std = np.array([ np.std(Q[:m]) for m in range(1,len(Q)) ])
        
        upperline = Avg + Std
        lowerline = Avg - Std

        x = np.arange(len(Avg))

        plt.plot(Avg,label="$\\left< Q \\right>$")
        plt.fill_between(x,lowerline,upperline,alpha=0.3,color='g',label="$\\sigma_Q$")
        plt.legend()
        plt.xlabel("Window size")
        plt.ylabel("$\\left< Q \\right>$")
        plt.title("Avg. Q over increasing window sizes")
        plt.savefig("avg_std_Q.pdf")
        plt.show()

if __name__ == "__main__":
    #T = "148.76"
    #Q = np.loadtxt(T+"_agg/Q.dat")
    #t, Econ, Etot = np.loadtxt("energyterms.xvg",usecols=(0,4,5),unpack=True)
    #t, Ebond, Eang = np.loadtxt("energyterms.xvg",usecols=(0,1,2),unpack=True)
    t,x = np.loadtxt("rmsd.xvg",unpack=True)
    plot_BSE(x)

