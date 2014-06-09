""" Plot loop lengths for all circular permutants


Description:

    The loop length of a contact represents the entropic cost of forming that
contact in the mean field approximation. The distribution of loop lengths 
therefore characterises the entropics of a particular contact map.


"""

import numpy as np
import matplotlib.pyplot as plt

def get_loops(M):
    loops = []
    for i in range(len(M)):
        for j in range(len(M)):
            if M[i,j] == 1:
                loops.append(abs(i-j))
    loops.sort()
    return np.array(loops)

def plot_loop_length_dist_vs_permutants():
    print "Loading Qref"
    Qref = np.loadtxt("Qref_cryst.dat")
    a = Qref
    Loops = []

    print "  Calculating loops of all permutants..."
    for k in range(len(Qref)):
        loops = get_loops(a)
        Loops.append(loops)
        a = np.roll(a,1,axis=0)
        a = np.roll(a,1,axis=1)

    Loops = np.array(Loops)
    bins = np.linspace(min(Loops.ravel()),max(Loops.ravel()),20)

    print "  Calculating loop length distributions..."
    Ldist = [ np.histogram(Loops[i,:],bins=bins)[0] for i in range(len(Loops)) ]
    Ldist = np.array(Ldist)

    avg_loops = [ np.mean(Loops[i,:]) for i in range(len(Loops)) ]
    std_loops = [ np.std(Loops[i,:]) for i in range(len(Loops)) ]

    ## Plot avg loop length and standard deviation
    #plt.plot(avg_loops,label="Avg. Loop")
    #plt.plot(std_loops,label="Std. Loop")
    #plt.legend()
    #plt.show()
    #raise SystemExit

    X,Y = np.meshgrid(bins,range(len(Qref)+1))

    print "Plotting"
    plt.pcolor(X,Y,Ldist)
    plt.xlim(min(bins),max(bins))
    plt.ylim(0,len(Qref))
    plt.xlabel("Loop Length")
    plt.ylabel("Permutant")
    cbar = plt.colorbar()
    cbar.set_label("Frequency")
    plt.title("Loop Length Distribution for all Permutants")
    print "Saving"
    plt.savefig("loop_length_dist_vs_permutant.pdf")
    plt.show()


if __name__ == "__main__":

    plot_loop_length_dist_vs_permutants()
