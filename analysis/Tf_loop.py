import numpy as np
import os

import contacts
import crunch_coordinates


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

        ## Calculates native contact reference matrix at lowest temperature.
        if (not os.path.exists(lowT+"/Qref_prob.dat")):
            os.chdir(lowT)
            Qref = contacts.probabilistic_reference()
            os.chdir(cwd2)
        else:
            Qref = np.loadtxt(lowT+"/Qref_prob.dat")

        ## Saves the reference matrix in each temp. directory. Submits PBS job for 
        ## calculating native contacts, native helical, and not-native contacts.
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            np.savetxt("Qref_prob.dat",Qref,delimiter=" ",fmt="%1d")
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


