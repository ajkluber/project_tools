import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt

def get_wham_config_basic(numbinsQ,minQ,stepQ,numbinsE,minE,stepE,temps):
    wham_config = "numDimensions 2 # number of dimensions in the dos histogram\n"
    wham_config += "               # energy and Q, energy is always in the first column, can\n"
    wham_config += "               # have up to 3 reaction coordinates\n\n"

    wham_config += "kB 0.008314    # Boltzmann constant\n\n"

    wham_config += "run_wham       # creates optimal density of states, writes to dosFile,\n"
    wham_config += "               # comment out to not run\n\n"

    wham_config += "dosFile dos    # density of states output filename for run_wham, input\n"
    wham_config += "               # filename for analysis (i.e. run_cv)\n\n"

    wham_config += "threads 1      # on multicore systems you can use up to 8 threads to speed\n"
    wham_config += "               # up the calculation\n\n"

    wham_config += "### energy binning ###\n"
    wham_config += "numBins %d\n" % numbinsE
    wham_config += "start %.2f\n" % minE
    wham_config += "step %.2f\n\n" % stepE

    wham_config += "### reaction coordinate 1 binning Q ###\n"
    wham_config += "numBins %d\n" % numbinsQ
    wham_config += "start %.2f\n" % minQ
    wham_config += "step %.2f\n\n" % stepQ

    wham_config += "### list of histogram filenames and their temperatures ###\n"
    wham_config += "numFiles %d\n" % len(temps)
    for temp in temps:
        wham_config += "name hist%d.dat temp %d\n" % temp
    wham_config += "\n"
    return wham_config

def get_wham_config_free_energy():
    wham_config += "#### Compute F(Q)\n"
    wham_config += "run_free            # comment out to not run, reads dosFile\n"
    wham_config += "run_free_out free   # prefix for the free energy curves\n"
    wham_config += "startTF 140         # first temperature to compute F(Q)\n"
    wham_config += "deltaTF 0.1         # step in temperature\n"
    wham_config += "ntempsF 300         # total F(Q) to generate\n\n"
    return wham_config

def get_wham_config_heat_capacity():
    wham_config += "#### Compute Cv(T)\n"
    wham_config += "run_cv              # comment out to not run, reads dosFile\n"
    wham_config += "run_cv_out cv       # filename for the temperature curves\n"
    wham_config += "startT 140          # starting temperature\n"
    wham_config += "deltaT 0.1          # step in temperature\n"
    wham_config += "ntemps 300          # total temps to generate\n\n"
    return wham_config

def get_wham_config_melting_curve():
    wham_config += "### Compute <Q>(T)\n"
    wham_config += "run_coord            # comment out to not run, reads dosFile\n"
    wham_config += "run_coord_out Q_vs_T # filename for the coordinate curve\n"
    wham_config += "startTC 140          # temperature to start at\n"
    wham_config += "deltaTC 0.1          # step in temperature\n"
    wham_config += "ntempsC 300          # total temps to generate\n\n"
    return wham_config
