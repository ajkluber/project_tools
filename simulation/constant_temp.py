""" Start simulations in the folding temperature loop. Tf_loop

Description:

    This module will be the library for submitting the simulation jobs for the
Tf_loop (folding temperature loop). The Tf_loop tries to find the folding
temperature by determining the melting curve over a large spread in
temperatures then narrowing in on the transition point. The goal is to obtain
equilibrium simulations at the folding temperature.

"""

import numpy as np
import subprocess as sb
from glob import glob
import os
import argparse
import shutil
import logging

import mdp

def main():
    """ Use gmxcheck on subdirectories.  """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--check', action='store_true', help='use gmxcheck on all subdirectories')
    args = parser.parse_args()
    if args.check == True:
        gmxcheck_subdirectories()
    else:
        pass

def check_completion(model,iteration,long=False):
    """ Checks to see if the previous Tf_loop simulation completed. 

    Description:

        First 
    checks the desired number of steps in the .mdp file then 
    checks to see if md.log has recorded that number of steps.
    """

    name = model.name

    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("simulation")

    cwd = os.getcwd()
    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir("%s/%s" % (cwd,sub))
    if long == True:
        tempfile = open("long_temps_last","r").readlines()
    else:
        tempfile = open("short_temps_last","r").readlines()
    #tempfile = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in tempfile  ]
    error = 0
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        os.chdir(tdir)
        print "  Running gmxcheck on ",tdir
        sb.call("gmxcheck -f traj.xtc",shell=True,stdout=open("check.out","w"),stderr=open("check.err","w")) 
        errorcode = "Fatal error"
        if (errorcode in open("check.err","r").read()) or (errorcode in open("check.out","r").read()):
            print "  FATAL ERROR in directory: ",tdir
            print "  somethings wrong with Gromacs traj.xtc file. See %s/check.err" % tdir
            error = 1

        # Determine the number of steps for completed run.
        for line in open("nvt.mdp" ,"r"):
            if line[:6] == "nsteps":
                nsteps = int(line.split()[2]) + 1
                break    
        finish_line = "Statistics over %d" % nsteps
        # Check if md.log has finished the required number of steps.
        if finish_line not in open("md.log","r").read():
            print "    Check %s simulation did not finish." % tdir
            # Restart run
            if os.path.exists("rst.pbs"):
                sb.call("qsub rst.pbs",shell=True)
                print "  Restarting: %s " % tdir
            error = 1
        os.chdir("..")

    if error == 1:
        print "  Cannot continue until simulations complete. Check if all unfinished runs were restarted properly."
        os.chdir(cwd)
        raise SystemExit
    else:
        if long == True:
            logger.info(" Finished: Equil_Tf")
        else:
            logger.info(" Finished: Tf_loop_iteration")
    os.chdir(cwd)

def gmxcheck_subdirectories():
    """ Run gmxcheck on all traj.xtc files in subdirecories. """
    runs = glob("*/traj.xtc")
    dirs = [ x[:-9] for x in runs ]
    error = 0
    for subdir in dirs:
        os.chdir(subdir)
        print "  Running gmxcheck on ",subdir
        sb.call("gmxcheck -f traj.xtc",shell=True,stdout=open("check.out","w"),stderr=open("check.err","w")) 
        errorcode = "Fatal error"
        if (errorcode in open("check.err","r").read()) or (errorcode in open("check.out","r").read()):
            print "  FATAL ERROR in directory: ",subdir
            print "  somethings wrong with Gromacs traj.xtc file. See %s/check.err" % subdir
            error = 1
        os.chdir("..")
    
    if error != 0:
        print "ERROR! Some trajectories did not pass gmxcheck."
        print " Exiting."
        raise SystemExit

def determine_new_temperatures():
    """ Find the temperatures which bracket the folding temperature.
        This takes the temperatures at which the average fraction of
        nonhelical contacts falls below 0.5 as bracketing the folding 
        temperature. A more complicated calculation is probably 
        needed for more complicated systems (proteins with intermediates)"""
    #temps = open("T_array_last.txt","r").readlines()
    temps = open("short_temps_last","r").readlines()
    temperatures = [ temp[:-1] for temp in temps ]
    temperatures.sort()
    lowerT = int(temperatures[0].split("_")[0])
    dT = int(temperatures[1].split("_")[0]) - lowerT
    upperT = int(temperatures[-1].split("_")[0])
    # Temperatures bracket the folding temperature if the average 
    # fraction of nonlocal native contacts goes from greater than 0.5 to less
    # than 0.5.
    for tdir in temperatures:
        Q = np.loadtxt(tdir+"/Qnhprob.dat")
        avgQ = np.mean(Q[len(Q)/2:])
        if avgQ > 0.5*max(Q):
            lowerT = int(tdir.split("_")[0])
        else:
            upperT = int(tdir.split("_")[0])
            break
        
    if dT == 2:
        # Previous run was final iteration. Now WHAM needs to be
        # done.
        newT_min, newT_max, newdeltaT = 0,0,0
    else:
        # Determine the finer grain T_array.
        newdeltaT = int(float(dT)/5.)
        # If newdeltaT < 5 then just do finest temperature spacing. 
        if newdeltaT < 5:
            newdeltaT = 2
            midT = int(0.5*(float(lowerT)+upperT))
            newT_min = midT - 20
            newT_max = midT + 20
        else:
            newT_min = lowerT + newdeltaT
            newT_max = upperT - newdeltaT
    #print "#DEBUGGING: New Ti, Tf, dT", newT_min, newT_max, newdeltaT
    return newT_min, newT_max, newdeltaT

def manually_extend_temperatures(model,iteration,method,temps,factor):
    """ To manually extend some temperatures """

    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("simulation")

    name = model.name
    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (model.path,name,iteration)
    # Determine directory to enter
    if method == "short":
        Tlist = [ str(int(x))+"_0" for x in temps ]
    elif method == "long":
        Tlist = []
        for i in range(len(temps)):
            for j in range(1,10):
                tempdir = "%.2f_%d" % (temps[i], j)
                if os.path.exists(sub+"/"+tempdir):
                    Tlist.append(tempdir) 
                else:
                    break

    os.chdir(sub)
    check_exist = [ os.path.exists(x) for x in Tlist ]
    # Check that all temperaures exist 
    if not all(check_exist):
        print "ERROR!"
        print "  Some temperature does not exist!"
        for i in range(len(Tlist)):
            print "  Temp:", Tlist[i], check_exist[i]  
        print "  Exiting."
        os.chdir(cwd)
        raise SystemExit

    cwd2 = os.getcwd()
    for k in range(len(Tlist)):
        Tdir = Tlist[k]
        os.chdir(Tdir)
        T = Tdir.split("_")[0]
        if model.dry_run == True:
            print "    Dryrun Success! " 
            os.chdir(cwd)
            raise SystemExit
        else:
            extend_temperature(T,factor)
            if os.path.exists("Q.dat"):
                os.remove("Q.dat")
            if os.path.exists("energyterms.xvg"):
                os.remove("energyterms.xvg")
        os.chdir(cwd2)
    if method == "short":
        logger.info(name,"Starting: Tf_loop_iteration")
    elif method == "long":
        logger.info(name,"Starting: Equil_Tf")

    os.chdir(cwd)

def extend_temperature(T,factor):
    """ Extend individual temperature run by factor """
    # Calculate new nsteps = factor*old_nsteps
    for line in open("nvt.mdp","r").readlines():
        if line.startswith("nsteps"):
            old_nsteps = int(line.split()[2])
            new_nsteps = str(int(round(factor*old_nsteps)))
            break
    
    # Save old .mdp and .tpr as something else.
    shutil.move("nvt.mdp","nvt.mdp")
    shutil.move("topol_4.6.tpr","old_topol_4.6.tpr")

    # Write new .mdp with more steps and recreate .tpr
    mdpfile = mdp.constant_temperature(T,new_nsteps)
    open("nvt.mdp","w").write(mdpfile)

    print "  Extending temp ", T, " to nsteps ",new_nsteps
    prep_step1 = 'grompp_sbm -n index.ndx -f nvt.mdp -c conf.gro -p topol.top -o topol_4.5.tpr '
    prep_step2 = 'grompp -n index.ndx -f nvt.mdp -c conf.gro -p topol.top -o topol_4.6.tpr '
    sb.call(prep_step1.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))
    sb.call(prep_step2.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))

    # Submit rst.pbs
    qsub = "qsub rst.pbs"
    sb.call(qsub.split(),stdout=open("rst.out","w"),stderr=open("rst.err","w"))

def folding_temperature_loop(model,iteration,new=False):
    """ The "folding temperature loop" is one of the several large-scale 
        logical structures in modelbuilder. It is entered anytime we want
        to determine the folding temperature. This could be when we have
        started a new project, refined the paramters, or returned to a 
        project in progress. The folding temperature loop successively 
        narrows in on the folding temperature."""

    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (model.path,model.name,iteration)
    if (not os.path.exists(sub)):
        os.mkdir(sub)
    # Check to see if the folding temperature has been found. If yes, then continue.
    if (not os.path.exists("%s/short_Tf" % sub)):
        os.chdir(sub)
        folding_temperature_loop_extension(model,new=new)
    else:
        # Folding temperature has been found. Continuing.
        pass
    os.chdir(cwd)

def folding_temperature_loop_extension(model,new=False):
    """ This is for doing an additional loop in the Tf_loop. It either starts
        an initial temperature array or refines the temperature range according
        to previous data. """
    # Check to see if a previous temperature range was used.
    #if (not os.path.exists("T_array_last.txt")) or new:
    if (not os.path.exists("short_temps_last.txt")) or new:
        # For initial exploration use very broad temperature increments.
        if model.initial_T_array != None:
            T_min = model.initial_T_array[0]
            T_max = model.initial_T_array[1]
            deltaT = model.initial_T_array[2]
        else:
            # Estimate folding temperature
            E = -model.native_stability
            N = float(model.n_residues)
            Tf_guess = int(round((36.081061*E/N) + 56.218196)) # calibration for LJ1210 contacts circa June 2014
            T_min = Tf_guess - 16
            T_max = Tf_guess + 16
            deltaT = 2
    else:
        # Use previous range to determine new range. 
        T_min, T_max, deltaT = determine_new_T_array()
    print "  Running temperature array: T_initial = %.2f   T_final = %.2f   dT = %.2f " % (T_min,T_max,deltaT)
    run_temperature_array(model,T_min,T_max,deltaT)
    logging.basicConfig(filename="%s.log" % model.name,level=logging.INFO)
    logger = logging.getLogger("simulation")
    logger.info("Submitting short_temps iteration ")
    logger.info("  T_min = %d , T_max = %d , dT = %d" % (T_min, T_max, deltaT))
    logger.info(" Starting: Tf_loop_iteration")

def start_next_Tf_loop_iteration(model,iteration):
    """ Estimate new folding temperature with calibration data

    Description:

        We made a calibration curve with the following points.
    """

    # Estimate folding temperature
    E = -model.native_stability
    N = float(model.n_residues)
    Tf_guess = (36.081061*E/N) + 56.218196 # calibration for LJ1210 contacts circa June 2014
    T_min = Tf_guess - 30
    T_max = Tf_guess + 30
    T_min = int(round(T_min))
    T_max = int(round(T_max))
    deltaT = 2
    name = model.name
    cwd = os.getcwd()
    sub = "%s/iteration_%d" % (name,iteration)
    if os.path.exists(sub):
        raise IOError("  The next iteration directory exists.")
    else:
        os.makedirs(sub)
    os.chdir(sub)
    run_temperature_array(model,T_min,T_max,deltaT)

    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("simulation")
    logger.info("Submitting T_array iteration %d" % iteration)
    logger.info("  T_min = %d , T_max = %d , dT = %d" % (T_min, T_max, deltaT))
    logger.info(" Starting: Tf_loop_iteration")
    os.chdir(cwd)

def manually_add_temperature_array(model,iteration,T_min,T_max,deltaT):
    """ To manually set the next temperature array."""
    logging.basicConfig(filename="%s.log" % model.name,level=logging.INFO)
    logger = logging.getLogger("simulation")
    sub = "%s/iteration_%d" % (model.name,iteration)
    os.chdir(sub)
    logger.info("Submitting T_array iteration %d " % iteration)
    logger.info("  T_min = %d , T_max = %d , dT = %d" % (T_min, T_max, deltaT))
    run_temperature_array(model,T_min,T_max,deltaT)
    logger.info(" Starting: Tf_loop_iteration")
    os.chdir("../..")

def manually_add_equilibrium_runs(model,iteration,temps):
    """ To manually set the next temperature array."""
    name = model.name
    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("simulation")

    sub = "%s/iteration_%d" % (name,iteration)
    os.chdir(sub)
    # Run for longer if the protein is really big.
    walltime, queue, ppn, nsteps = determine_equil_walltime(model)

    T_string = ''
    for i in range(len(temps)):
        T = "%.2f" % temps[i]
        for simnum in range(1,4):
            simpath = T+"_"+str(simnum)
            # Only start the simulation if directory doesn't exist.
            if (not os.path.exists(simpath)):
                T_string += "%s\n" % simpath
                os.mkdir(simpath)
                os.chdir(simpath)
                print "    Running temperature ", simpath
                run_constant_temp(model,T,name,nsteps=nsteps,walltime=walltime,queue=queue)
                os.chdir("..")
            else:
                # Directory exists for this temperature: continue.
                continue

    open("long_temps","a").write(T_string)
    open("long_temps_last","w").write(T_string)
    logger.info(" Starting: Equil_Tf")
    os.chdir("../..")

def run_equilibrium_simulations(model,iteration):
    """ Run very long (equilibrium) simulations at the estimated folding 
        temperature."""

    name = model.name
    cwd = os.getcwd()
    sub = "%s/%s/iteration_%d" % (cwd,name,iteration)
    Tf = open("%s/short_Tf" % sub ,"r").read().split()[0]

    walltime, queue, ppn, nsteps = determine_equil_walltime(model)

    os.chdir(sub)
    T_string = ''
    for n in range(3):
        T = "%.2f" % (float(Tf)+1.*(n-1))
        for simnum in range(1,4):
            simpath = T+"_"+str(simnum)
            # Only start the simulation if directory doesn't exist.
            if (not os.path.exists(simpath)):
                T_string += "%s\n" % simpath
                os.mkdir(simpath)
                os.chdir(simpath)
                print "    Running temperature ", simpath
                run_constant_temp(model,T,name,nsteps=nsteps,walltime=walltime,queue=queue)
                os.chdir("..")
            else:
                # Directory exists for this temperature: continue.
                continue

    open("long_temps","a").write(T_string)
    open("long_temps_last","w").write(T_string)
    logging.basicConfig(filename="%s.log" % name,level=logging.INFO)
    logger = logging.getLogger("simulation")
    logger.info(" Starting: Equil_Tf")
    os.chdir(cwd)

def determine_equil_walltime(model):
    """ Estimate an efficient walltime."""
    N = model.n_residues
    ppn = "1"
    nsteps = "500000000"
    queue="serial"
    if N < 60:
        walltime="16:00:00"
    else:
        if N > 100:
            queue="serial_long"
            walltime="20:00:00"
        else:
            walltime="15:00:00"
    return walltime, queue, ppn,nsteps

def determine_walltime(model):
    """ Estimate an efficient walltime."""
    N = model.n_residues
    ppn = "1"
    nsteps = "100000000"
    queue="serial"
    if N < 60:
        walltime="03:00:00"
    else:
        if N > 100:
            if N > 200:
                walltime="16:00:00"
            else:
                walltime="5:00:00"
        else:
            walltime="3:00:00"
    return walltime, queue, ppn,nsteps

def run_temperature_array(model,T_min,T_max,deltaT):
    """ Simulate range of temperatures to find the folding temperature. """

    Temperatures = range(T_min,T_max+deltaT,deltaT)
    # Run for longer if the protein is really big.
    walltime, queue, ppn, nsteps = determine_walltime(model)
    name = model.name
    T_string = ''
    for T in Temperatures:
        simpath = str(T)+"_0"
        # Only start the simulation is directory doesn't exist.
        if (not os.path.exists(simpath)):
            T_string += "%d_0\n" % T
            os.mkdir(simpath)
            os.chdir(simpath)
            if not model.dry_run:
                print "  Running temperature ", T
            run_constant_temp(model,T,name,nsteps=nsteps,walltime=walltime,queue=queue,ppn=ppn)
            os.chdir("..")
        else:
            continue
    open("short_temps","a").write(T_string)
    open("short_temps_last","w").write(T_string)
    open("short_Ti_Tf_dT.txt","w").write("%d %d %d" % (T_min, T_max, deltaT))

def run_constant_temp(model,T,name,nsteps="100000000",walltime="23:00:00",queue="serial",ppn="1"):
    """ Start a constant temperature simulation with Gromacs. 

    Description:

        Save the grofile and topology file in model object, then 
    submit a pbs job to run simulation.

    """
    # Loading and writing grompp.
    mdpfile = mdp.constant_temperature(str(T),nsteps)
    open("nvt.mdp","w").write(mdpfile)

    # Write all needed simulation files.
    model.save_simulation_files()

    # Start simulation
    jobname = "%s_%s" % (name,str(T))
    prep_run(jobname,walltime=walltime,queue=queue,ppn=ppn)

    if model.dry_run == True:
        print "    Dryrun: Successfully saved simulation files for ", T
    else:
        qsub = "qsub run.pbs"
        sb.call(qsub.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))
    
def get_pbs_string(jobname,queue,ppn,walltime):
    """ Return basic PBS job script. """
    pbs_string = "#!/bin/bash \n"
    pbs_string +="#PBS -N %s \n" % jobname
    pbs_string +="#PBS -q %s \n" % queue
    pbs_string +="#PBS -l nodes=1:ppn=%s \n" % ppn
    pbs_string +="#PBS -l walltime=%s \n" % walltime
    pbs_string +="#PBS -V \n\n"
    pbs_string +="cd $PBS_O_WORKDIR\n"
    pbs_string +="mdrun -s topol_4.6.tpr"
    return pbs_string

def get_rst_pbs_string(jobname,queue,ppn,walltime):
    """ Return basic PBS job script for restarting. """
    pbs_string = "#!/bin/bash \n"
    rst_string = "#!/bin/bash \n"
    rst_string +="#PBS -N %s_rst \n" % jobname
    rst_string +="#PBS -q %s \n" % queue
    rst_string +="#PBS -l nodes=1:ppn=%s \n" % ppn
    rst_string +="#PBS -l walltime=%s \n" % walltime
    rst_string +="#PBS -V \n\n"
    rst_string +="cd $PBS_O_WORKDIR\n"
    rst_string +="mdrun -s topol_4.6.tpr -cpi state.cpt"
    return rst_string

def prep_run(jobname,walltime="23:00:00",queue="serial",ppn="1"):
    """ Executes the constant temperature runs."""

    prep_step1 = 'grompp_sbm -n index.ndx -f nvt.mdp -c conf.gro -p topol.top -o topol_4.5.tpr'
    prep_step2 = 'grompp -n index.ndx -f nvt.mdp -c conf.gro -p topol.top -o topol_4.6.tpr'
    sb.call(prep_step1.split(),stdout=open("prep.out","w"),stderr=open("prep.err","w"))
    sb.call(prep_step2.split(),stdout=open("prep.out","w"),stderr=open("prep.err","w"))

    pbs_string = get_pbs_string(jobname,queue,ppn,walltime)
    open("run.pbs","w").write(pbs_string)

    rst_string = get_rst_pbs_string(jobname,queue,ppn,walltime)
    open("rst.pbs","w").write(rst_string)

if __name__ == '__main__':
    main()
