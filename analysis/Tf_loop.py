import numpy as np
import os

import contacts
import crunch_coordinates
import wham


def analyze_temperature_array(System,i,append_log):
    ''' Analyze the previously simulated temperatures of a Tf_loop iteration.
        Goes into the active Tf_loop directory and crunches all coordinates.
        Exits after submitting a couple PBS scripts to compute Q and 
        energyterms.xvg.
    '''
    cwd = os.getcwd()
    if System.error[i] == 0:
        sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
        os.chdir(cwd+"/"+sub)
        tempfile = open("T_array_last.txt","r").readlines()
        temperatures = [ temp[:-1] for temp in tempfile  ]
        allTs = [ temp[:-1] for temp in open("T_array.txt","r").readlines() ]
        allTs.sort()
        lowT = allTs[0]

        cwd2 = os.getcwd()
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            if (not os.path.exists("rmsd.xvg")) or (not os.path.exists("radius_cropped.xvg")) or \
               (not os.path.exists("energyterms.xvg")) or (not os.path.exists("phis.xvg")):
                crunch_coordinates.crunch_all(System.subdirs[i]+"_"+tdir)
            os.chdir(cwd2)

        ## Calculates native contact reference matrix at lowest temperature. NOT USED FOR GO MODEL.
        if (not os.path.exists(lowT+"/Qref_cryst.dat")):
            print "## DEBUGGING: Calculating reference matrix."
            print "## DEBUGGING: ",lowT
            os.chdir(lowT)
            Qref = contacts.probabilistic_reference()
            os.chdir(cwd2)
        else:
            Qref = np.loadtxt(lowT+"/Qref_cryst.dat")

        ## Saves the reference matrix in each temp. directory. Submits PBS job for 
        ## calculating native contacts, native helical, and not-native contacts.
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            np.savetxt("Qref_cryst.dat",Qref,delimiter=" ",fmt="%1d")
            crunch_coordinates.crunch_Q(System.subdirs[i]+"_"+tdir)

            os.chdir(cwd2)
        os.chdir(cwd)
        append_log(System.subdirs[i],"Starting: Tf_loop_analysis")
        System.append_log(System.subdirs[i],"Starting: Tf_loop_analysis")
    else:
        pass

def check_completion(System,i,append_log):
    ''' Check if the Tf_loop_analysis finished by seeing if all needed files
        were generated.
    '''
    done = 0
    cwd = os.getcwd()
    sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    tempfile = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in tempfile  ]
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir(cwd2+"/"+tdir)
        if  os.path.exists("rmsd.xvg") and os.path.exists("radius_cropped.xvg") and \
            os.path.exists("energyterms.xvg") and os.path.exists("phis.xvg") and \
            os.path.exists("Qprob.dat"):
            System.append_log(System.subdirs[i],"  analysis done for "+tdir)
            done = 1
        else:
            done = 0
        os.chdir(cwd2)

    os.chdir(cwd)
    if done == 1:
        append_log(System.subdirs[i],"Finished: Tf_loop_analysis")
        System.append_log(System.subdirs[i],"Finished: Tf_loop_analysis")

def check_if_wham_is_next(System,i,append_log):
    ''' Check if the last temperature step, dT=1. If it was start 
        prepping and running WHAM calculation for the Heat Capacity.'''

    cwd = os.getcwd()
    sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    Tinfo = open("Ti_Tf_dT.txt","r").read().split()
    Ti,Tf,dT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])

    if dT == 2:
        ## Its time for WHAM
        print "Temperature interval has reached dT=2. Time for WHAM."
        print "Starting wham_Cv..."
        System.append_log(System.subdirs[i],"  prepping wham_Cv inputs")
        if os.path.exists(cwd2+"/wham"):
            pass
        else:
            os.makedirs(cwd2+"/wham")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"HeatCap")

        System.append_log(System.subdirs[i],"  running wham for heat capacity")
        append_log(System.subdirs[i],"Starting: wham_Cv")
        wham.run_wham("HeatCap")
        flag = 1
    else:
        ## Its not time for WHAM
        flag = 0

    os.chdir(cwd)
    return flag

def continue_wham(System,i,append_log):
    ''' If WHAM has already run for the Heat Capacity (Cv) then prep files
        and run WHAM for 1D & 2D free energy surfaces.'''

    cwd = os.getcwd()
    sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    Tinfo = open("Ti_Tf_dT.txt","r").read().split()
    Ti,Tf,dT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])

    ## Check for completion
    if os.path.exists(cwd2+"/wham/Heat_rmsd_Rg.dat"):
        print "Finished wham_Cv..."
        System.append_log(System.subdirs[i],"  wham heat capacity done")
        append_log(System.subdirs[i],"Finished: wham_Cv")
        print "Starting wham_FreeEnergy..."
        append_log(System.subdirs[i],"Starting: wham_FreeEnergy")
        System.append_log(System.subdirs[i],"  prepping wham inputs for 1D PMFs")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"1DFreeEnergy")
        System.append_log(System.subdirs[i],"  running wham for 1D PMFs")
        wham.run_wham("1DFreeEnergy")
        System.append_log(System.subdirs[i],"  prepping wham inputs for 2D PMFs")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"FreeEnergy")
        System.append_log(System.subdirs[i],"  running wham for 2D PMFs")
        wham.run_wham("FreeEnergy")
    else:
        print "wham_Cv may not have finished. Check if Heat_rmsd_Rg.dat exists."
        print "Exiting."

    os.chdir(cwd)
