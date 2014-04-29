import numpy as np
import matplotlib.pyplot as plt

def plot_kinetic_mechanism():
    ''' The kinetic mechanism is defined as the mechanism of only successful or
        "reactive" trajectories (i.e. those that traverse from unfolded to
        folded before returning to unfolded). This trims away unreactive 
        fluctuactions from each of the states.

        NOT DONE.'''
    coord = "Q"
    statefile = open(coord+"_states.txt","r").readlines()[1:]
    bounds = []
    for line in statefile:
        bounds.append([line.split()[0],float(line.split()[1]),float(line.split()[2]),float(line.split()[3])])

    print "Loading Q.dat, Qres.dat"
    Q = np.loadtxt("Q.dat")
    Qres = np.loadtxt("Qres.dat")


    minQ = min(Q)
    maxQ = max(Q)
    bins = 50
    incQ = (float(maxQ) - float(minQ))/bins

    Qprogress = np.zeros((bins,len(Qres[0])),float)
    counts = np.zeros(bins,float)

    left_bound = bounds[0][3]
    right_bound = bounds[2][2]
    folding = 0
    for i in range(1,len(Q)):

        if folding == 0:
            if (Q[i] > left_bound) and (Q[i-1] < left_bound):
                folding = 1

def plot_thermodynamic_mechanism():

    print "Loading Q.dat, Qres.dat"
    Q = np.loadtxt("Q.dat")
    Qres = np.loadtxt("Qres.dat")

    minQ = min(Q)
    maxQ = max(Q)
    bins = 50
    incQ = (float(maxQ) - float(minQ))/bins

    Qprogress = np.zeros((bins,len(Qres[0])),float)
    counts = np.zeros(bins,float)


    print "Histogram Qres depending on Q"
    for i in range(len(Q)):
        for n in range(bins):
            if ((minQ + n*incQ)  < Q[i]) and (Q[i] < (minQ + (n+1)*incQ)):
                Qprogress[n,:] += Qres[i,:]
                counts[n] += 1

    native = Qres[0,:]
    native[ native == 0 ] = 1

    print "Calculating fraction of formed native contacts for each bin"
    Qprogress = ((Qprogress/native).T)/counts

    print "Plotting thermodynamic mechanism"
    plt.figure()
    plt.pcolor(Qprogress,edgecolors='k')
    cbar = plt.colorbar()
    cbar.set_label("Fraction local contacts formed $Q_{local}$")
    plt.xlabel("Folding Progress from $[Q_{min},Q_{max}] = [%.2f,%.2f]$" % (minQ/float(maxQ),1.0))
    plt.ylabel("Sequence index")
    plt.title("Thermodynamic Folding Progress")
    plt.savefig("mechanism_profile.pdf")
    plt.show()


if __name__ == "__main__":

    plot_thermodynamic_mechanism()





