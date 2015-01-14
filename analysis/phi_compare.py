import numpy as np
import matplotlib.pyplot as plt

from phi_ss import *


def get_exp_phi_ss(name):
    element, bounds = get_sec_structure(name)

    lines = [ x.rstrip("\n") for x in open("%s/mutants/core.ddG" % name,"r").readlines() ]

    exp_phi_ss = [ [] for i in range(len(element)) ]

    for line in lines[1:]:
        vals = line.split()
        if int(vals[8]) == 1:
            continue
        else:
            mut_idx = int(vals[0])
            for j in range(len(element)):
                if (mut_idx >= bounds[j][0]) and (mut_idx <= bounds[j][1]):
                    exp_phi_ss[j].append(float(vals[7]))
                    break

    #print exp_phi_ss    ## DEBUGGIN
    avg_phi_ss = np.zeros(len(element))
    for j in range(len(element)):
        avg_phi_ss[j] = np.mean(exp_phi_ss[j])

    return avg_phi_ss,element

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of protein to plot.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to plot.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration

    sim_phi_ss = calculate_average_phi_ss(name,iteration)

    exp_phi_ss,element = get_exp_phi_ss(name) 

    #print sim_phi_ss
    #print exp_phi_ss

    ## Calculate least squares fit and R^2 value.
    x = np.array(sim_phi_ss)
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A,exp_phi_ss)[0]
    fit = m*x + c*np.ones(len(x))

    SEy = np.sum((exp_phi_ss - np.mean(exp_phi_ss))**2)
    SEline = np.sum((exp_phi_ss - fit)**2)
    r2 = 1. - (SEline/SEy)
     
    ## Plot
    plt.plot([0,1],[0,1],'k',lw=2)
    #plt.plot(x,fit,'b',lw=2,label="y = %.3fx + %.3f  R$^2$=%.3f" % (m,c,r2))
    plt.plot(x,fit,'b',lw=2,label="R$^2$=%.3f" % r2)
    for i in range(len(element)):
        plt.plot(sim_phi_ss[i],exp_phi_ss[i],'o',label=element[i])
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.grid(True)
    plt.legend()
    plt.xlabel("simulation $\phi_{ss}$")
    plt.ylabel("experiment $\phi_{ss}$")
    plt.title("$\phi_{ss}$ comparison for %s iteration %d" % (name,iteration))
    if not os.path.exists("%s/iteration_%d/plots" % (name,iteration)):
        os.mkdir("%s/iteration_%d/plots" % (name,iteration))
    plt.savefig("%s/iteration_%d/plots/compare_phi_ss.png" % (name,iteration))
    plt.savefig("%s/iteration_%d/plots/compare_phi_ss.pdf" % (name,iteration))
    plt.show()
