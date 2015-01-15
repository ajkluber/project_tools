''' Run Weighted Histogram Analysis Method Smog tool.

Description:

    Functions to prepare input files and run WHAM to compute 
thermal averages. For details on WHAM see reference (1). This
code uses a WHAM executable that is part of the SMOG set of 
tools, see reference (2) for more information .



(1) Kumar,S.; Bouzida,D.; Swendsen, R.H.; Kollman, P.A.; Rosenberg, J.
The Weighted Histogram Analysis Method for Free-energy calculations on
biomolecules. I. The Method. J. Comput. Chem. 1992, 13, 1011-1021.

(2) Noel, J.K.; Whitford, P.C.; Sanbonmatsu, K.Y.; Onuchic, J.N.
SMOG@ctbp: Simplified Deployment of Structure-Based Models in GROMACS.
Nucleic Acids Res. 2010, 38, W657-61.

'''



import matplotlib.pyplot as plt
import numpy as np
import subprocess as sb
import os
import shutil

def get_wham_config_basic(numbinsE,minE,stepE,numbinsC,minC,stepC,dim):
    wham_config = "numDimensions %d # number of dimensions in the dos histogram\n" % dim
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
    wham_config += "start %.5f\n" % minE
    wham_config += "step %.5f\n\n" % stepE

    wham_config += "### reaction coordinate 1 binning ###\n"
    wham_config += "numBins %d\n" % numbinsC
    wham_config += "start %.5f\n" % minC
    wham_config += "step %.5f\n\n" % stepC

    return wham_config

def get_wham_config_coordinate_binning(numbinsC,minC,stepC):
    wham_config = "### additional reaction coordinate binning ###\n"
    wham_config += "numBins %d\n" % numbinsC
    wham_config += "start %.5f\n" % minC
    wham_config += "step %.5f\n\n" % stepC
    return wham_config

def get_wham_config_coordinate_output(outfile,T):
    wham_config = "### Output reaction coordinate binning ###\n"
    wham_config += "run_coord \n"
    wham_config += "run_coord_out %s\n" % outfile
    wham_config += "startTC %.5f\n\n" % T
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

def run_wham_expdH_k(mut,Tf,bounds):
    ''' Prepare histogram files for wham.
    
        Concatenates all the data from the same temperature for the histogram
    files in the wham subdirectory.
    '''
    ## To Do:
    ##  1. Run WHAM to get F(Q) and Cv(T). DONE
    ##      1a. User adjusts T for F(Q) to find Tf --> Save as Mut_0/whamQ/Tf.txt
    ##      1b. User determines bounds for TS, U, N states --> Save as Mut_0/whamQ/state_bounds
    ##  2. Calculate dH  TESTING --> WORKS
    ##  3. Run WHAM to get the following thermal averages as function
    ##    of Q at Tf:
    ##    3a. <Q_ij>_X the interaction energy of contact ij.
    ##    3b. <exp(-beta*dH_k)>_X the perturbation due to mutation k.
    ##    [ 3c. <Q_i>_X the contact probability ]
    ##
    ##  4. Calculate DeltaDeltaG's for simulation and phi_values
    ##  5. Calculate matrix M and solve for new parameters.

    directories = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ]
    temperatures = [ x.split("_")[0] for x in open("T_array_last.txt","r").readlines() ]
    start_at = 1

    unique_temps = []
    counts = []
    for t in temperatures:
        if t not in unique_temps:
            unique_temps.append(t)
            counts.append(temperatures.count(t))
        else:
            pass

    ## Concatenate the E and Q data at a particular temperature then
    ## save in the whamQ subdirectory
    cwd = os.getcwd()
    histfilenames = ""
    histfiles = []
    i = 0  
    print "  Preparing histograms for ", mut
    for k in range(len(unique_temps)):
        T = unique_temps[k]
        beta = 1./float(T)
        histname = "dH_hist_%s_%.2f" % (mut,float(T))
        #if not os.path.exists("whamQ/"+histname):
        if True:
            print "    temperature: ",T
            for n in range(1,counts[k]+1):
                
                q = np.loadtxt(T+"_"+str(n)+"/Q.dat")
                state_indicator = np.zeros(len(q),int)
                ## Assign every frame a state label. State indicator is integer 1-N for N states.
                for state_num in range(len(bounds)-1):
                    instate = (q > bounds[state_num]).astype(int)*(q <= bounds[state_num+1]).astype(int)
                    state_indicator[instate == 1] = state_num+1
                if any(state_indicator == 0):
                    print "ERROR! There are unassigned frames!"
                    print sum((state_indicator == 0).astype(int)), " unassigned frames out of ", len(q)
                    print " Exiting"
                    raise SystemExit
                
                expdH = np.exp(-beta*np.loadtxt(T+"_"+str(n)+"/dH_"+mut+".dat"))
                eng = np.loadtxt(T+"_"+str(n)+"/energyterms.xvg",usecols=(4,5))[:,1]
                if n == 1:
                    State_indicator = state_indicator
                    ExpdH = expdH
                    E = eng
                else:
                    State_indicator = np.concatenate((State_indicator,state_indicator),axis=0)
                    ExpdH = np.concatenate((ExpdH,expdH),axis=0)
                    E = np.concatenate((E,eng),axis=0)
            if i == 0:
                maxexpdH = max(ExpdH)
                minexpdH = min(ExpdH)
                maxE = max(E)
                minE = min(E[1:])
            else:
                maxexpdH = max([max(ExpdH),maxexpdH])
                minexpdH = min([min(ExpdH),minexpdH])
                maxE = max([max(E),maxE])
                minE = min([min(E[1:]),minE])
            histogram = np.zeros((len(ExpdH),3),float)
            histogram[:,0] = E
            histogram[:,1] = State_indicator
            histogram[:,2] = ExpdH
            if long:
                np.savetxt("long_wham/"+histname,histogram,fmt="%15.6f")
            else:
                np.savetxt("short_wham/"+histname,histogram,fmt="%15.6f")
        else:
            print " loading ",histname
            if long:
                E,state,ExpdH = np.loadtxt("long_wham/"+histname,unpack=True)
            else:
                E,state,ExpdH = np.loadtxt("short_wham/"+histname,unpack=True)
            if i == 0:
                maxexpdH = max(ExpdH)
                minexpdH = min(ExpdH)
                maxE = max(E)
                minE = min(E[1:])
            else:
                maxexpdH = max([max(ExpdH),maxexpdH])
                minexpdH = min([min(ExpdH),minexpdH])
                maxE = max([max(E),maxE])
                minE = min([min(E[1:]),minE])
        if long:
            histfiles.append("long_wham/"+histname)
        else:
            histfiles.append("short_wham/"+histname)
        histfilenames += "name %s temp %.2f\n" % (histname,float(T))
        i += 1

    print "  histogram files: ",histfiles
    if long:
        os.chdir("long_wham")
    else:
        os.chdir("short_wham")
    ## Binning settings
    numbinsState = len(bounds)-1
    startState = 0.2
    stepState = 1.0

    stepE = 5 
    numbinsE = int(round((maxE - minE)/stepE))

    startExpdH = 0.9
    numbinsExpdH = 50
    stepExpdH = (maxexpdH - startExpdH)/float(numbinsExpdH)

    wham_basic = get_wham_config_basic(numbinsE,minE,stepE,numbinsState,startState,stepState,3)
    wham_basic += get_wham_config_coordinate_binning(numbinsExpdH,startExpdH,stepExpdH)
    wham_basic += get_wham_config_coordinate_output("expdH_"+mut,Tf)
    wham_basic += "### list of histogram filenames and their temperatures ###\n"
    wham_basic += "numFiles %d\n" % len(unique_temps)
    wham_basic += histfilenames
    wham_basic += "\n"

    #print wham_basic       ## DEBUGGING

    open("%s.config" % mut,"w").write(wham_basic)

    print "  Running wham for ", mut
    cmd1 = "java -jar WHAM.1.06.jar --config %s.config" % mut
    sb.call(cmd1.split(),stdout=open(mut+".out","w"),stderr=open(mut+".err","w"))

    os.chdir("..")
    #return wham_basic,temperatures

def prepare_histograms_heat_capacity(long=False):
    ''' Prepare histogram files for wham.
    
        Concatenates all the data from the same temperature for the histogram
    files in the wham subdirectory.
    '''
    if long:
        if not os.path.exists("long_wham"):
            os.mkdir("long_wham")
        directories = [ x.rstrip("\n") for x in open("long_temps_last","r").readlines() ]
        temperatures = [ x.split("_")[0] for x in open("long_temps_last","r").readlines() ]
        endings = [ int(x.split("_")[1].rstrip("\n")) for x in open("long_temps_last","r").readlines() ]
    else:
        if not os.path.exists("short_wham"):
            os.mkdir("short_wham")
        directories = [ x.rstrip("\n") for x in open("short_temps_last","r").readlines() ]
        temperatures = [ x.split("_")[0] for x in open("short_temps_last","r").readlines() ]
        endings = [ int(x.split("_")[1].rstrip("\n")) for x in open("short_temps_last","r").readlines() ]

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
    ## save in the wham subdirectory
    cwd = os.getcwd()
    histfilenames = ""
    i = 0  
    print "  Preparing histograms..."
    for k in range(len(unique_temps)):
        T = unique_temps[k]

        print "    temperature: ",T
        for n in range(start_at,counts[k]+1):
            q = np.loadtxt("%s_%d/Q.dat" % (T,n))
            eng = np.loadtxt("%s_%d/energyterms.xvg" % (T,n),usecols=(5,))
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
        if long:
            np.savetxt("long_wham/hist_%.2f" % float(T),histogram,fmt="%15.6f")
        else:
            np.savetxt("short_wham/hist_%.2f" % float(T),histogram,fmt="%15.6f")
        histfilenames += "name hist_%.2f temp %.2f\n" % (float(T),float(T))
        i += 1

    if long == True:
        stepQ = 2 
        stepE = 2 
    else:
        stepQ = 5 
        stepE = 5 
    numbinsQ = int(round((maxQ - minQ)/stepQ))
    numbinsE = int(round((maxE - minE)/stepE))
    wham_basic = get_wham_config_basic(numbinsE,minE,stepE,numbinsQ,minQ,stepQ,2)
    wham_basic += "### list of histogram filenames and their temperatures ###\n"
    wham_basic += "numFiles %d\n" % len(unique_temps)
    wham_basic += histfilenames
    wham_basic += "\n"

    return wham_basic,temperatures

def run_wham_for_heat_capacity(model,long=False):
    ''' Prepare wham histograms and run Jeff's WHAM code'''

    wham_basic,temperatures = prepare_histograms_heat_capacity(long=long)

    temps = [ float(x) for x in temperatures ]
    startT = min(temps)
    if long == True:
        startT -= 2.
        stopT = max(temps) + 2.
        deltaTCv = 0.01
        ntempsCv = int(round((stopT - startT)/deltaTCv))
        ntempsF = 30
        #ntempsF = (float(max(temps)) - float(startT))/deltaTF
        deltaTF = (float(max(temps)) - float(startT))/ntempsF
    else:
        stopT = max(temps)
        deltaTCv = 0.05
        deltaTF = 0.4
        ntempsCv = int(round((stopT - startT)/deltaTCv))
        ntempsF = int(round((stopT - startT)/deltaTF))

    wham_Cv = get_wham_config_heat_capacity(startT,deltaTCv,ntempsCv)
    wham_Melt = get_wham_config_melting_curve(startT,deltaTCv,ntempsCv)
    wham_Free = get_wham_config_free_energy(startT,deltaTF,ntempsF)
    
    if long:
        os.chdir("long_wham")
    else:
        os.chdir("short_wham")
    #PROJECTS = os.environ["PROJECTS"]
    PROJECT_TOOLS = os.environ["PROJECT_TOOLS"]
    shutil.copy("%s/project_tools/analysis/WHAM.1.06.jar" % PROJECT_TOOLS,".")

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

    ax2.plot(QvsT[:,0],QvsT[:,1]/model.n_pairs,'b')
    ax2.set_ylim(0,1)
    ax2.set_ylabel("$\\left< Q \\right>(T)$")
    plt.savefig("cv_and_melt.pdf")
    Tf = Cv[list(Cv[:,1]).index(max(Cv[:,1])),0]
    print "  Folding temperature: ", Tf
    os.chdir("..")
    if long == True:
        print "  Wham done! Plotted Cv and melting curve: %s/iteration_%d/long_wham/cv_and_melt.pdf" % (model.subdir,model.iteration)
        open("long_Tf","w").write("%.2f" % Tf)
    else:
        print "  Wham done! Plotted Cv and melting curve: %s/iteration_%d/short_wham/cv_and_melt.pdf" % (model.subdir,model.iteration)
        open("short_Tf","w").write("%.2f" % Tf)
