import numpy as np
import subprocess as sb
import os
import shutil

import contacts
import crunch_coordinates
import wham
import plot


def aggregate_equilibrium_runs(System,append_log,reagg=False):
    ''' Aggregate equilibrium simulation data into one directory for    
        ease of access.'''

    append_log(System.subdir,"Starting: Aggregating_Equil_Runs")
    cwd = os.getcwd()
    sub = System.subdir+"/"+System.mutation_active_directory
    os.chdir(cwd+"/"+sub)
    temps = [ x.split('_')[0] for x in open("T_array.txt","r").readlines() ] 
    unique_temps = []
    counts = []
    for t in temps:
        if t not in unique_temps:
            unique_temps.append(t)
            counts.append(temps.count(t))
        else:
            pass

    coords = ["Q.dat","A.dat","Qh.dat","Qnh.dat","Qlocal.dat","Qnonlocal.dat","Nh.dat",
              "Qres.dat","Qhres.dat","Qnhres.dat","Qlocalres.dat","Qnonlocalres.dat",
              "rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg",
              "Qprob.dat","Qhprob.dat","Qnhprob.dat"]

    for k in range(len(unique_temps)):
        T = unique_temps[k]
        if not os.path.exists(T+"_agg"):
            os.mkdir(T+"_agg")
        print "  Aggregating data into directory %s_agg" % T

        if (not os.path.exists(T+"_agg/traj.xtc")) or reagg:
            print "  Concatenating trajectories for ", T
            xtcfiles = ''
            for n in range(1,counts[k]+1):
                xtcfiles += " "+T+"_"+str(n)+"/traj.xtc "
            cmd1 = "trjcat -f "+xtcfiles+" -o "+T+"_agg/traj.xtc -cat"
            sb.call(cmd1,shell=True,stderr=open(T+"_agg/trjcat.err","w"),stdout=open(T+"_agg/trjcat.out","w"))

        shutil.copy(T+"_1/Native.pdb",T+"_agg/Native.pdb")
        shutil.copy(T+"_1/BeadBead.dat",T+"_agg/BeadBead.dat")

        for cord in coords:
            if (not os.path.exists(T+"_agg/"+cord)) or reagg:
                print "    Aggregating", cord
                for i in range(1,counts[k]+1):
                    x = np.loadtxt(T+"_"+str(i)+"/"+cord)
                    if i == 1:
                        X = x
                    else:
                        X = np.concatenate((X,x),axis=0)
                np.savetxt(T+"_agg/"+cord, X)
         
    os.chdir(cwd)
    append_log(System.subdir,"Finished: Aggregating_Equil_Runs")

def analyze_temperature_array(System,append_log,equil=False):
    ''' Analyze the previously simulated temperatures of a Tf_loop iteration.
        Goes into the active Tf_loop directory and crunches all coordinates.
        Exits after submitting a couple PBS scripts to compute Q and 
        energyterms.xvg.
    '''
    cwd = os.getcwd()
    if System.error == 0:
        if equil == True:
            sub = System.subdir+"/"+System.mutation_active_directory
            qwalltime = "00:08:00"
            cwalltime = "00:10:00"
        else:
            sub = System.subdir+"/"+System.Tf_active_directory
            qwalltime = "00:04:00"
            cwalltime = "00:02:00"
        print "  Analyzing temperature array for", sub
        os.chdir(cwd+"/"+sub)
        tempfile = open("T_array_last.txt","r").readlines()
        temperatures = [ temp[:-1] for temp in tempfile  ]
        allTs = [ temp[:-1] for temp in open("T_array.txt","r").readlines() ]
        allTs.sort()
        lowT = allTs[0]

        cwd2 = os.getcwd()
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            os.chdir(cwd2+"/"+tdir)
            crunchfiles = ["rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg",
                            "Qprob.dat","Qhprob.dat","Qnhprob.dat"]
            flag = all([ os.path.exists(file) for file in crunchfiles ])
            if not flag:
                print "    Crunching coordinates for ",tdir
                crunch_coordinates.crunch_all(System.subdir+"_"+tdir,walltime=cwalltime)
                crunch_coordinates.crunch_Q(System.subdir+"_"+tdir,walltime=qwalltime)
            else:
                print "    Skipping directory ",tdir
            os.chdir(cwd2)

        ## Saves the reference matrix in each temp. directory. Submits PBS job for 
        ## calculating native contacts, native helical, and not-native contacts.
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            os.chdir(cwd2)
        os.chdir(cwd)
        if equil == True:
            append_log(System.subdir,"Starting: Equil_Tf_analysis")
            System.append_log("Starting: Equil_Tf_analysis")
        else:
            append_log(System.subdir,"Starting: Tf_loop_analysis")
            System.append_log("Starting: Tf_loop_analysis")
    else:
        pass

def check_completion(System,append_log,equil=False):
    ''' Check if the Tf_loop_analysis finished by seeing if all needed files
        were generated.
    '''
    done = 0
    cwd = os.getcwd()
    if equil == True:
        sub = System.subdir+"/"+System.mutation_active_directory
        qwalltime = "00:20:00"
        cwalltime = "00:10:00"
    else:
        sub = System.subdir+"/"+System.Tf_active_directory
        qwalltime = "00:04:00"
        cwalltime = "00:02:00"
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    print "  Checking analysis in directory "+sub
    tempfile = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in tempfile  ]
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir(cwd2+"/"+tdir)
        if  os.path.exists("rmsd.xvg") and os.path.exists("radius_cropped.xvg") and \
            os.path.exists("energyterms.xvg") and os.path.exists("phis.xvg") and \
            os.path.exists("Qprob.dat"):
            if os.path.exists("Nh.dat"): 
                print "    Crunch coordinates done. "
            else:
                print "    Crunch coordinates done. Crunching Nh for "+tdir
                System.append_log("    crunching Nh for "+tdir)
                crunch_coordinates.crunch_Nh()
            System.append_log("    analysis done for "+tdir)
            done = 1
        else:
            print "    Crunching not done. Retrying for "+tdir
            crunch_coordinates.crunch_all(System.subdir+"_"+tdir,walltime=cwalltime)
            crunch_coordinates.crunch_Q(System.subdir+"_"+tdir,walltime=qwalltime)
            done = 0
        os.chdir(cwd2)

    os.chdir(cwd)
    if done == 1:
        print "  Analysis completed."
        if equil == True:
            append_log(System.subdir,"Finished: Equil_Tf_analysis")
            System.append_log("Finished: Equil_Tf_analysis")
        else:
            append_log(System.subdir,"Finished: Tf_loop_analysis")
            System.append_log("Finished: Tf_loop_analysis")
    else:
        print "  Analysis has not finished."

def check_if_wham_is_next(System,append_log):
    ''' Check if the last temperature step, dT=1. If it was start 
        prepping and running WHAM calculation for the Heat Capacity.'''

    cwd = os.getcwd()
    sub = System.subdir+"/"+System.Tf_active_directory
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    Tinfo = open("Ti_Tf_dT.txt","r").read().split()
    Ti,Tf,dT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])

    if dT == 2:
        ## Its time for WHAM
        print "Temperature interval has reached dT=2. Time for WHAM."
        print "Starting wham_Cv..."
        System.append_log("  prepping wham_Cv inputs")
        if os.path.exists(cwd2+"/wham"):
            pass
        else:
            os.makedirs(cwd2+"/wham")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"HeatCap")

        System.append_log("  running wham for heat capacity")
        append_log(System.subdir,"Starting: wham_Cv")
        wham.run_wham("HeatCap")
        flag = 1
    else:
        ## Its not time for WHAM
        flag = 0

    os.chdir(cwd)
    return flag

def continue_wham(System,append_log):
    ''' If WHAM has already run for the Heat Capacity (Cv) then prep files
        and run WHAM for 1D & 2D free energy surfaces.'''

    cwd = os.getcwd()
    sub = System.subdir+"/"+System.Tf_active_directory
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    Tinfo = open("Ti_Tf_dT.txt","r").read().split()
    Ti,Tf,dT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])

    ## Check for completion
    if os.path.exists(cwd2+"/wham/Heat_rmsd_Rg.dat"):
        print "Finished wham_Cv..."
        System.append_log("  wham heat capacity done")
        append_log(System.subdir,"Finished: wham_Cv")
        print "Starting wham_FreeEnergy..."
        append_log(System.subdir,"Starting: wham_FreeEnergy")
        System.append_log("  prepping wham inputs for 1D PMFs")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"1DFreeEnergy")
        System.append_log("  running wham for 1D PMFs")
        wham.run_wham("1DFreeEnergy")
        System.append_log("  prepping wham inputs for 2D PMFs")
        wham.prep_input_files(Ti,Tf,dT,cwd2,"FreeEnergy")
        System.append_log("  running wham for 2D PMFs")
        wham.run_wham("FreeEnergy")
    else:
        print "wham_Cv may not have finished. Check if Heat_rmsd_Rg.dat exists."
        print "Exiting."

    os.chdir(cwd)
