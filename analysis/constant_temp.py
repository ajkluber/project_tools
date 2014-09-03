import numpy as np
import subprocess as sb
import os
import shutil

import contacts
import crunch_coordinates
import wham
import plot

def aggregate_equilibrium_runs(model,append_log,reagg=False):
    ''' Aggregate equilibrium simulation data into one directory for    
        ease of access.

    DEPRECATED 7-6-2014
    '''

    append_log(model.subdir,"Starting: Aggregating_Equil_Runs")
    cwd = os.getcwd()
    sub = model.subdir+"/"+model.mutation_active_directory
    os.chdir(cwd+"/"+sub)
    temperatures = [ x.split('_')[0] for x in open("T_array.txt","r").readlines() ] 
    unique_temps = []
    counts = []
    for t in temperatures:
        if t not in unique_temps:
            unique_temps.append(t)
            counts.append(temperatures.count(t))
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
    append_log(model.subdir,"Finished: Aggregating_Equil_Runs")

def determine_walltime(model,long=False):
    ''' Estimate an efficient walltime.'''
    N = model.n_residues
    ppn = "4"
    queue = "serial"
    if N < 60:
        qwalltime = "00:04:00"
        cwalltime = "00:02:00"
    else:
        if N > 160:
            queue = "bigmem"
            if N > 250:
                if long:
                    qwalltime = "00:16:00"
                    cwalltime = "00:10:00"
                else:
                    qwalltime = "00:10:00"
                    cwalltime = "00:08:00"
                ppn = "8"
            else:
                if long:
                    qwalltime = "00:16:00"
                    cwalltime = "00:10:00"
                else:
                    qwalltime = "00:12:00"
                    cwalltime = "00:10:00"
                ppn = "6"
        else:
            qwalltime = "00:08:00"
            cwalltime = "00:10:00"
    return qwalltime, cwalltime, ppn, queue

def analyze_temperature_array(model,append_log,equil=False):
    ''' Analyze the previously simulated temperatures of a Tf_loop iteration.
        Goes into the active Tf_loop directory and crunches all coordinates.
        Exits after submitting a couple PBS scripts to compute Q and 
        energyterms.xvg.
    '''
    cwd = os.getcwd()
    if model.error == 0:
        if equil == True:
            sub = model.subdir+"/Mut_"+str(model.Mut_iteration)
            qwalltime = "00:10:00"
            cwalltime = "00:07:00"
        else:
            sub = model.subdir+"/Tf_"+str(model.Tf_iteration)
            qwalltime = "00:05:00"
            cwalltime = "00:03:00"
        ppn = "1"
        queue = "serial"
        print "  Analyzing temperatures in", sub
        os.chdir(cwd+"/"+sub)
        temperatures = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ]

        cwd2 = os.getcwd()
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            os.chdir(cwd2+"/"+tdir)
            crunchfiles = ["rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg"]
            crunchQfiles = ["Q.dat","qimap.dat"]

            flag = all([ os.path.exists(file) for file in crunchfiles ])
            flagQ = all([ os.path.exists(file) for file in crunchQfiles ])
            if (not flag) or (not flag):
                if not flag:
                    print "    Crunching coordinates for ",tdir
                    crunch_coordinates.crunch_all(model.subdir+"_"+tdir,walltime=cwalltime,ppn=ppn)
                if not flagQ:
                    print "    Crunching Q for ",tdir
                    crunch_coordinates.crunch_Q(model.subdir+"_"+tdir,walltime=qwalltime,ppn=ppn,queue=queue)
            else:
                print "    Skipping directory ",tdir
            os.chdir(cwd2)

        os.chdir(cwd)
        if equil == True:
            append_log(model.subdir,"Starting: Equil_Tf_analysis")
            append_log("Starting: Equil_Tf_analysis",subdir=True)
        else:
            append_log(model.subdir,"Starting: Tf_loop_analysis")
            append_log("Starting: Tf_loop_analysis",subdir=True)
    else:
        pass

def check_completion(model,append_log,equil=False):
    ''' Check if the Tf_loop_analysis finished by seeing if all needed files
        were generated.
    '''
    done = 0
    cwd = os.getcwd()
    if equil == True:
        sub = model.subdir+"/Mut_"+str(model.Mut_iteration)
        qwalltime = "00:10:00"
        cwalltime = "00:06:00"
    else:
        sub = model.subdir+"/Tf_"+str(model.Tf_iteration)
        qwalltime = "00:05:00"
        cwalltime = "00:03:00"
    os.chdir(cwd+"/"+sub)
    cwd2 = os.getcwd()
    print "  Checking analysis in directory "+sub
    temperatures = [ x.rstrip("\n") for x in open("T_array_last.txt","r").readlines() ]
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir(cwd2+"/"+tdir)
        files = ["rmsd.xvg","radius_gyration.xvg","energyterms.xvg","phis.xvg","Q.dat","qimap.dat"] 
        check_files = all([ os.path.exists(file) for file in files ])
        if check_files:
            #if os.path.exists("Nh.dat"): 
            #    print "    Crunch coordinates done. "
            #else:
            #    print "    Crunch coordinates done. Crunching Nh for "+tdir
            #    append_log("    crunching Nh for "+tdir,subdir=True)
            #    crunch_coordinates.crunch_Nh()
            print "    Saving Qh, Qnh, Qlocal, Qnonlocal for %s" % tdir
            append_log("    Saving Qh, Qnh, Qlocal, Qnonlocal for %s" % tdir,subdir=True)
            crunch_coordinates.reorganize_qimap()
            append_log("    analysis done for %s" % tdir,subdir=True)
            done = 1
        else:
            print "    Crunching not done. Retrying for "+tdir
            crunch_coordinates.crunch_all("%s_%s" % (model.subdir,tdir),walltime=cwalltime)
            crunch_coordinates.crunch_Q("%s_%s" % (model.subdir,tdir),walltime=qwalltime)
            done = 0
        os.chdir(cwd2)

    os.chdir(cwd)
    if done == 1:
        print "  Analysis completed."
        if equil == True:
            append_log(model.subdir,"Finished: Equil_Tf_analysis")
            append_log("Finished: Equil_Tf_analysis",subdir=True)
        else:
            append_log(model.subdir,"Finished: Tf_loop_analysis")
            append_log("Finished: Tf_loop_analysis",subdir=True)
    else:
        print "  Analysis has not finished."

def run_wham_heat_capacity(model,append_log,Mut=False):
    ''' Check if the last temperature step, dT=1. If it was start 
        prepping and running WHAM calculation for the Heat Capacity.'''

    cwd = os.getcwd()
    print "*** NOTE: module load jdk required for WHAM ***"
    if Mut == True:
        sub = model.subdir+"/Mut_"+str(model.Mut_iteration)
        os.chdir(cwd+"/"+sub)
        cwd2 = os.getcwd()
        print "Running wham for heat capacity, free energy curves, and melting curve"
        if not os.path.exists("whamQ"):
            os.mkdir("whamQ")
        append_log("  running wham for heat capacity, free energy, and melting curve",subdir=True)
        append_log(model.subdir,"Starting: Equil_Tf_wham")
        wham.run_wham_for_heat_capacity(model,Mut=True)
        append_log(model.subdir,"Finished: Equil_Tf_wham")
        flag = 1
    else:
        sub = model.subdir+"/Tf_"+str(model.Tf_iteration)
        os.chdir(cwd+"/"+sub)
        cwd2 = os.getcwd()
        Tinfo = open("Ti_Tf_dT.txt","r").read().split()
        T_min,T_max,deltaT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])
        if deltaT <= 5:
            ## Its time for WHAM
            print "Since deltaT <=5 --> Time for WHAM."
            print "Running wham for heat capacity, free energy curves, and melting curve"
            if not os.path.exists("whamQ"):
                os.mkdir("whamQ")
            append_log("  running wham for heat capacity, free energy, and melting curve",subdir=True)
            append_log(model.subdir,"Starting: Tf_wham")
            wham.run_wham_for_heat_capacity(model)
            append_log(model.subdir,"Finished: Tf_wham")
            flag = 1
        else:
            ## Its not time for WHAM
            flag = 0

    os.chdir(cwd)
    return flag

