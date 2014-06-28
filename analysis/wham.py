import matplotlib.pyplot as plt
import numpy as np
import subprocess as sb
import os
import shutil

def get_wham_config_basic(numbinsQ,minQ,stepQ,numbinsE,minE,stepE):
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
    return wham_config

def get_wham_config_heat_capacity(startT,deltaT,ntemps):
    wham_config = "#### Compute Cv(T)\n"
    wham_config += "run_cv              # comment out to not run, reads dosFile\n"
    wham_config += "run_cv_out cv       # filename for the temperature curves\n"
    wham_config += "startT %6.2f       # starting temperature\n"     % startT
    wham_config += "deltaT %6.2f       # step in temperature\n"      % deltaT
    wham_config += "ntemps %6d        # total temps to generate\n\n" % ntemps
    return wham_config

def get_wham_config_free_energy(startTF,deltaTF,ntempsF):
    wham_config = "#### Compute F(Q)\n"
    wham_config += "run_free            # comment out to not run, reads dosFile\n"
    wham_config += "run_free_out free   # prefix for the free energy curves\n"
    wham_config += "startTF %6.2f      # first temperature to compute F(Q)\n" % startTF
    wham_config += "deltaTF %6.2f      # step in temperature\n"           % deltaTF
    wham_config += "ntempsF %6d      # total F(Q) to generate\n\n"          % ntempsF
    return wham_config

def get_wham_config_melting_curve(startTC,deltaTC,ntempsC):
    wham_config = "### Compute <Q>(T)\n"
    wham_config += "run_coord            # comment out to not run, reads dosFile\n"
    wham_config += "run_coord_out Q_vs_T # filename for the coordinate curve\n"
    wham_config += "startTC %6.2f       # temperature to start at\n"   % startTC
    wham_config += "deltaTC %6.2f       # step in temperature\n"       % deltaTC
    wham_config += "ntempsC %6d        # total temps to generate\n\n" % ntempsC
    return wham_config

def prepare_histograms(Mut=False):
    """ Prepare histogram files for wham.
    
        Concatenates all the data from the same temperature for the histogram
    files in the whamQ subdirectory.
    """
    if not os.path.exists("whamQ"):
        os.mkdir("whamQ")

    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ]
    temperatures = [ x.split("_")[0] for x in open("T_array_last.txt","r").readlines() ]
    endings = [ int(x.split("_")[1].rstrip("\n")) for x in open("T_array_last.txt","r").readlines() ]
    start_at = min(endings)

    unique_temps = []
    counts = []
    for t in temperatures:
        if t not in unique_temps:
            unique_temps.append(t)
            if start_at == 0:
                counts.append(temperatures.count(t)-1)
            else:
                counts.append(temperatures.count(t))
        else:
            pass

    ## Concatenate the E and Q data at a particular temperature then
    ## save in the whamQ subdirectory
    cwd = os.getcwd()
    histfilenames = ""
    i = 0  
    print "  Preparing histograms..."
    for k in range(len(unique_temps)):
        T = unique_temps[k]

        print "    temperature: ",T
        for n in range(start_at,counts[k]+1):
            q = np.loadtxt(T+"_"+str(n)+"/Q.dat")
            eng = np.loadtxt(T+"_"+str(n)+"/energyterms.xvg",usecols=(4,5))[:,1]
            if n == start_at:
                Q = q
                E = eng
            else:
                Q = np.concatenate((Q,q),axis=0)
                E = np.concatenate((E,eng),axis=0)
        if i == 0:
            maxQ = max(Q)
            minQ = min(Q)
            maxE = max(E)
            minE = min(E)
        else:
            maxQ = max([max(Q),maxQ])
            minQ = min([min(Q),minQ])
            maxE = max([max(E),maxE])
            minE = min([min(E),minE])

        histogram = np.zeros((len(Q),2),float)
        histogram[:,0] = E
        histogram[:,1] = Q
        np.savetxt("whamQ/hist_%.2f" % float(T),histogram,fmt="%15.6f")
        histfilenames += "name hist_%.2f temp %.2f\n" % (float(T),float(T))
        i += 1

    if Mut == True:
        stepQ = 2 
        stepE = 2 
    else:
        stepQ = 5 
        stepE = 5 
    numbinsQ = int(round((maxQ - minQ)/stepQ))
    numbinsE = int(round((maxE - minE)/stepE))
    wham_basic = get_wham_config_basic(numbinsQ,minQ,stepQ,numbinsE,minE,stepE)
    wham_basic += "numFiles %d\n" % len(unique_temps)
    wham_basic += histfilenames
    wham_basic += "\n"

    return wham_basic,temperatures

def run_wham_for_heat_capacity(model,Mut=False):
    """ Prepare wham histograms and run Jeff's WHAM code"""

    wham_basic,temperatures = prepare_histograms(Mut=Mut)

    temps = [ float(x) for x in temperatures ]
    startT = min(temps)
    if Mut == True:
        deltaTCv = 0.01
        ntempsCv = int(round((float(max(temps)) - float(startT))/deltaTCv))
        ntempsF = 10
        deltaTF = (float(max(temps)) - float(startT))/ntempsF
    else:
        deltaTCv = 0.05
        deltaTF = 0.4
        ntempsCv = int(round((float(max(temps)) - float(startT))/deltaTCv))
        ntempsF = int(round((float(max(temps)) - float(startT))/deltaTCv))

    wham_Cv = get_wham_config_heat_capacity(startT,deltaTCv,ntempsCv)
    wham_Melt = get_wham_config_melting_curve(startT,deltaTCv,ntempsCv)
    wham_Free = get_wham_config_free_energy(startT,deltaTF,ntempsF)
    
    os.chdir("whamQ")
    PROJECTS = os.environ["PROJECTS"]
    shutil.copy(PROJECTS+"/project_tools/analysis/WHAM.1.06.jar",".")

    open("cv.config","w").write(wham_basic + wham_Cv)
    open("free.config","w").write(wham_basic + wham_Free)
    open("melt.config","w").write(wham_basic + wham_Melt)

    cmd1 = "java -jar WHAM.1.06.jar --config cv.config"
    cmd2 = "java -jar WHAM.1.06.jar --config free.config"
    cmd3 = "java -jar WHAM.1.06.jar --config melt.config"

    print "  Running wham..."
    sb.call(cmd1.split(),stdout=open("cv.out","w"),stderr=open("cv.err","w"))
    sb.call(cmd2.split(),stdout=open("free.out","w"),stderr=open("free.err","w"))
    sb.call(cmd3.split(),stdout=open("melt.out","w"),stderr=open("melt.err","w"))

    Cv = np.loadtxt("cv",usecols=(0,1))
    QvsT = np.loadtxt("Q_vs_T",dtype=float)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(Cv[:,0],Cv[:,1],'r')
    ax1.set_xlabel("Temperature (K)")
    ax1.set_ylabel("Heat Capacity (kJ/mol K)")
    ax1.set_title("$C_v(T)$ and $\\left< Q \\right>(T)$ for %s" % model.name)

    ax2.plot(QvsT[:,0],QvsT[:,1]/model.n_contacts,'b')
    ax2.set_ylim(0,1)
    ax2.set_ylabel("$\\left< Q \\right>(T)$")
    plt.savefig("cv_and_melt.pdf")
    if Mut == True:
        print "  Wham done! Plotted Cv and melting curve: %s/Mut_%d/whamQ/cv_and_melt.pdf" % (model.subdir,model.Mut_iteration)
    else:
        print "  Wham done! Plotted Cv and melting curve: %s/Tf_%d/whamQ/cv_and_melt.pdf" % (model.subdir,model.Tf_iteration)
    Tf = Cv[list(Cv[:,1]).index(max(Cv[:,1])),0]
    print "  Folding temperature: ", Tf
    os.chdir("..")
    open("Tf.txt","w").write("%.2f" % Tf)
