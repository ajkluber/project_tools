import numpy as np
import subprocess as sb
import os
import shutil

import crunch_coordinates
import wham

def determine_walltime(model,long=False):
    """ Estimate an efficient walltime."""
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

def analyze_temperature_array(model,iteration,long=False):
    """ Analyze the previously simulated temperatures of a Tf_loop iteration.
        Goes into the active Tf_loop directory and crunches all coordinates.
        Exits after submitting a couple PBS scripts to compute Q and 
        energyterms.xvg.
    """
    name = model.name
    cwd = os.getcwd()
    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir("%s/%s" % (cwd,sub))
    if long:
        temperatures = [ x.rstrip("\n") for x in open("long_temps","r").readlines() ]
        qwalltime = "00:01:00"
        cwalltime = "00:01:00"
    else:
        temperatures = [ x.rstrip("\n") for x in open("short_temps","r").readlines() ]
        qwalltime = "00:01:00"
        cwalltime = "00:01:00"
    ppn = "1"
    queue = "serial"
    print "  Analyzing temperatures in", sub

    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir("%s/%s/%s" % (cwd,sub,tdir))
        #crunchfiles = ["rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg"]
        #crunchQfiles = ["Q.dat","qimap.dat"]
        crunchfiles = ["energyterms.xvg"]
        crunchQfiles = ["Q.dat"]

        flag = all([ os.path.exists(file) for file in crunchfiles ])
        flagQ = all([ os.path.exists(file) for file in crunchQfiles ])
        if (not flag) or (not flag):
            print "    Calculating energies and Q for ",tdir
            crunch_coordinates.crunch_all("%s_%s" % (name,tdir),model.contact_type,walltime=cwalltime,ppn=ppn,n_tables=model.n_tables)
            crunch_coordinates.crunch_Q("%s_%s" % (name,tdir),model.contact_type,walltime=qwalltime,ppn=ppn,queue=queue)
        else:
            print "    Skipping directory ",tdir
        os.chdir("%s/%s" % (cwd,sub))

    os.chdir(cwd)
    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("analysis")
    if long == True:
        logger.info(" Starting: Equil_Tf_analysis")
        logger.info(" Starting: Equil_Tf_analysis",subdir=True)
    else:
        logger.info(" Starting: Tf_loop_analysis")
        logger.info(" Starting: Tf_loop_analysis",subdir=True)

def check_completion(model,iteration,long=False):
    """ Check if the Tf_loop_analysis finished by seeing if all needed files
        were generated.
    """
    name = model.name
    done = 0
    cwd = os.getcwd()
    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir("%s/%s" % (cwd,sub))
    if long == True:
        temperatures = [ x.rstrip("\n") for x in open("long_temps","r").readlines() ]
        qwalltime = "00:01:00"
        cwalltime = "00:01:00"
    else:
        temperatures = [ x.rstrip("\n") for x in open("short_temps","r").readlines() ]
        qwalltime = "00:01:00"
        cwalltime = "00:01:00"
    print "  Checking analysis in directory %s" % sub
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir("%s/%s/%s" % (cwd,sub,tdir))
        #files = ["rmsd.xvg","radius_gyration.xvg","energyterms.xvg","phis.xvg","Q.dat","qimap.dat"] 
        files = ["Q.dat","energyterms.xvg"] 
        check_files = all([ os.path.exists(file) for file in files ])
        if check_files:
            done = 1
        else:
            print "    Crunching not done. Retrying for %s" % tdir
            crunchfiles = ["rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg"]
            crunchQfiles = ["Q.dat","qimap.dat"]

            flag = all([ os.path.exists(file) for file in crunchfiles ])
            flagQ = all([ os.path.exists(file) for file in crunchQfiles ])

            if not flag:
                crunch_coordinates.crunch_all("%s_%s" % (model.name,tdir),model.contact_type,walltime=cwalltime,n_tables=model.n_tables)

            if not flagQ:
                crunch_coordinates.crunch_Q("%s_%s" % (model.name,tdir),model.contact_type,walltime=qwalltime)
            done = 0
        os.chdir("%s/%s" % (cwd,sub))

    os.chdir(cwd)
    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("analysis")
    if done == 1:
        print "  Analysis completed."
        if long == True:
            logger.info(" Finished: Equil_Tf_analysis")
        else:
            logger.info(" Finished: Tf_loop_analysis")
    else:
        print "  Analysis has not finished."
        raise SystemExit

def run_wham_heat_capacity(model,iteration,long=False):
    """ Check if the last temperature step, dT=1. If it was start 
        prepping and running WHAM calculation for the Heat Capacity."""

    name = model.name
    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("analysis")
    cwd = os.getcwd()
    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir("%s/%s" % (cwd,sub))
    print "*** NOTE: module load jdk required for WHAM ***"
    if long:
        print "Running wham for heat capacity, free energy curves, and melting curve"
        if not os.path.exists("long_wham"):
            os.mkdir("long_wham")
        wham.run_wham_for_heat_capacity(long=True)
        logger.info(" Starting: Equil_Tf_wham")
        logger.info(" Finished: Equil_Tf_wham")
        flag = 1
    else:
        Tinfo = open("short_Ti_Tf_dT.txt","r").read().split()
        T_min,T_max,deltaT = int(Tinfo[0]), int(Tinfo[1]), int(Tinfo[2])
        if deltaT <= 5:
            ## Its time for WHAM
            print "Since deltaT <=5 --> Time for WHAM."
            print "Running wham for heat capacity, free energy curves, and melting curve"
            if not os.path.exists("short_wham"):
                os.mkdir("short_wham")
            wham.run_wham_for_heat_capacity()
            logger.info(" Starting: Tf_wham")
            logger.info(" Finished: Tf_wham")
            flag = 1
        else:
            ## Its not time for WHAM
            flag = 0

    os.chdir(cwd)
    return flag

