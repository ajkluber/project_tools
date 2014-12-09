""" Calculate Transition Path Theory (TPT) quantities.


    Partitions the trajectory into transition paths (TP's) and not TP's. TP's
are stretches of the trajectory that are "reactive" in that they involve a
successful transition between two states. 

    Transition Path Theory (TPT) then has some useful formulas for calculating 
useful dynamical quantities of interest, like rates.

    One quantity computed is the probability of being on a transition path
given the value of a reaction coordiante (r). We can compute this easily using
Bayes' Theorem as follows:

                    P(TP|r) = P(r|TP)*P(TP)/P(r)        eq. (9) in ref (1)

    Where,
    P(TP|r) = Probability of being on TP given value of coordinate r.
    P(r|TP) = Probability of having value r on a TP.
    P(TP) = Probability of being on a TP in entire trajectory.
    P(r)  = Probability of having value r in entire trajectory.

To Do:
    - Compute more interesting things

References:
(1) Hummer; "From Transition paths to transition states and rate coefficients".
    J. Chem. Phys. 2004. vol 120. number 2

"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_state_bounds():
    """ Bounds for each state. Bounds are bin edges along Q. """
    if os.path.exists("state_bounds.txt"):
        statefile = open("state_bounds.txt","r").readlines()
    else:
        print "ERROR!"
        print "  Please create state_bounds.txt"
        print "  With the boundaries of each state along Q"
        print "  Exiting"
        raise SystemExit
    
    state_bounds = []
    state_labels = []
    for line in statefile:
        info = line.split()
        state_labels.append(info[0])
        state_bounds.append(float(info[1]))
        state_bounds.append(float(info[2]))
    
    return state_bounds,state_labels

def partition_TP(Q,stateA,stateB):
    """ Partition the trajectory into transition paths (TP) and not TP."""

    
    if (Q[0] < stateA):
        prevState = "A"
        indicator = [0]
    elif (Q[0] > stateB):
        prevState = "B"
        indicator = [2]
    else:
        prevState = "TP"
        indicator = [1]

    TP = []
    notTP = []
    tempTP = []
    TPlengths = []
    for i in range(1,len(Q)):

        if (Q[i] < stateA):
            currState = "A" 
            temp = 0
        elif (Q[i] > stateB):
            currState = "B" 
            temp = 2
        else:
            currState = "TP" 
            temp = 1

        indicator.append(temp)

        if (currState == "TP") and ((prevState == "A") or (prevState == "B")):
            ## If you are starting a potential TP
            tempTP = [Q[i]]
            fromState = prevState
        elif ((currState == "A") or (currState == "B")) and (prevState == "TP"):
            ## If you are ending a potential TP
            if currState == fromState:
                notTP.extend(tempTP)
            else:
                ## Reactive trajectory!
                #print fromState, currState, prevState, len(tempTP)     # DEBUGGING
                TP.extend(tempTP)
                TPlengths.append(float(len(tempTP)))
            fromState = currState
        elif ((currState == "A") and (prevState == "A")) or \
             ((currState == "B") and (prevState == "B")):
            notTP.append(Q[i])
        elif ((currState == "TP") and (prevState == "TP")):
            tempTP.append(Q[i])
        else:
            print "never"

        prevState = currState

    indicator = np.array(indicator)
    n_TPs = len(TPlengths)

    return TP,notTP,TPlengths,n_TPs,indicator
        
def TP_analysis(TP,notTP,TPlengths,indicator):
    """ Analyze the transition paths from simulation

    """

    n_TPs = len(TPlengths)

    hA = float((indicator == 0).astype(int))       # Indicator function for state A
    hB = float((indicator == 2).astype(int))       # Indicator function for state B
    
    totalframes = hA + hB
    cA = hA/totalframes                 # Fraction of state A
    cB = hB/totalframes                 # Fraction of state B

    tau_TP = np.mean(TPlengths)          # Average transition path length

    maxQ = float(max(Q))
    minQ = float(min(Q))
    bins = np.linspace(minQ,maxQ,n_bins)

    N = float(len(Q))                           # Num frames total
    N_TP = float(len(TP))                       # Num frames on TP's
    Nr,bins = np.histogram(Q,bins=bins)         # Num times reaction coord falls in bin
    Nr_TP,bins = np.histogram(TP,bins=bins)     # Num times reaction coord falls in bin on TP

    P_TP = N_TP/N                       # Prob of being on TP
    Pr = Nr.astype(float)/N             # Prob of having reaction coord value
    Pr_TP = Nr_TP.astype(float)/N_TP    # Prob of having reaction coord value on TP
    P_TP_r = Pr_TP*P_TP/Pr              # Prob of being on TP given reaction coord value

    k1 = P_TP/(2.*cA*tau_TP)            # Folding rate
    k2 = P_TP/(2.*cB*tau_TP)            # Unfolding rate

def transition_state_diffusion_constant():
    """ Calculate the diffusion coefficient at the transition state.
        
    """
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--name', type=str, required=True, help='Name of subdirectory.')
    parser.add_argument('--iteration', type=int, required=True, help='Iteration to calculate.')
    parser.add_argument('--directory', type=str, required=True, help='Directory.')
    parser.add_argument('--coordinate', type=str, default="Q.dat", help='Name of coordinate file to use.')
    parser.add_argument('--bins', type=int, default=30, help='Num of bins along coordinate.')
    args = parser.parse_args()

    name = args.name
    iteration = args.iteration
    directory = args.directory
    coord = args.coordinate
    n_bins = args.bins

    cwd = os.getcwd()

    os.chdir("%s/iteration_%d" % (name,iteration))
    state_bounds, state_labels = get_state_bounds()
    stateA = state_bounds[1]
    stateB = state_bounds[4]

    if not os.path.exists(directory):
        print "ERROR! Directory %s does not exist!" % directory
        print " Exiting"
        raise SysteExit

    T = directory.split("_")[0]
    dir = directory

    ## Patition trajectory into reactive transition paths based on state
    ## boundaries along some coordinate Q.
    os.chdir(dir)
    Q = np.loadtxt(coord)
    TP, notTP, TPlengths, n_TPs, indicator = partition_TP(Q,stateA,stateB)

    ## Do TP analysis and save interesting data.
    if not os.path.exists("TPT"):
        os.mkdir("TPT")
    os.chdir("TPT")
    #TP_analysis(TP,notTP,TPlengths,indicator)
    n_TPs = len(TPlengths)

    hA = float(sum((indicator == 0).astype(int)))       # Indicator function for state A
    hB = float(sum((indicator == 2).astype(int)))       # Indicator function for state B
    
    totalframes = hA + hB
    cA = hA/totalframes                 # Fraction of state A
    cB = hB/totalframes                 # Fraction of state B

    tau_TP = np.mean(TPlengths)          # Average transition path length

    maxQ = float(max(Q))
    minQ = float(min(Q))
    bins = np.linspace(minQ,maxQ,n_bins)

    N = float(len(Q))                           # Num frames total
    N_TP = float(len(TP))                       # Num frames on TP's
    Nr,bins = np.histogram(Q,bins=bins)         # Num times reaction coord falls in bin
    Nr_TP,bins = np.histogram(TP,bins=bins)     # Num times reaction coord falls in bin on TP

    P_TP = N_TP/N                       # Prob of being on TP
    Pr = Nr.astype(float)/N             # Prob of having reaction coord value
    Pr_TP = Nr_TP.astype(float)/N_TP    # Prob of having reaction coord value on TP
    P_TP_r = Pr_TP*P_TP/Pr              # Prob of being on TP given reaction coord value

    k1 = P_TP/(2.*cA*tau_TP)            # Folding rate
    k2 = P_TP/(2.*cB*tau_TP)            # Unfolding rate

    os.chdir("..")
    os.chdir("..")

    os.chdir(cwd)
