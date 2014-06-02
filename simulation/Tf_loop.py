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

def check_completion(System,append_log,equil=False):
    """ Checks to see if the previous Tf_loop simulation completed. First 
        checks the desired number of steps in the grompp.mdp file then 
        checks to see if md.log has recorded that number of steps."""
    cwd = os.getcwd()
    if equil == True:
        sub = System.subdir+"/"+System.mutation_active_directory
    else:
        sub = System.subdir+"/"+System.Tf_active_directory
    os.chdir(cwd+"/"+sub)
    tempfile = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in tempfile  ]
    error = 0
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        check_error = run_gmxcheck(tdir)
        ## Determine the number of steps for completed run.
        for line in open(tdir+"/"+"grompp.mdp","r"):
            if line[:6] == "nsteps":
                nsteps = int(line.split()[2]) + 1
                break    
        finish_line = "Statistics over " + str(nsteps)
        ## Check if md.log has finished the required number of steps.
        if finish_line in open(tdir+"/md.log","r").read():
            System.append_log("  %s finished." % tdir)
        else:
            print "    Check %s simulation did not finish." % tdir
            #print "    Cannot continue with errors."
            ## Try to restart the run if possible.
            if os.path.exists(tdir+"/rst.pbs"):
                os.chdir(tdir)
                qrst = "qsub rst.pbs"
                sb.call(qrst.split(),stdout=open("rst.out","w"),stderr=open("rst.err","w"))
                os.chdir(cwd+"/"+sub)
                System.append_log("  %s did not finish. restarting" % tdir)
                print "  %s did not finish. Restarting: submitting rst.pbs. " % tdir
            else:
                System.append_log("  %s did not finish. did not find a rst.pbs. skipping." % tdir)
            error = 1

    if error == 1:
        print "  Cannot continue until simulations complete. Check if all unfinished runs were restarted properly."
        pass 
    else:
        if equil == True:
            append_log(System.subdir,"Finished: Equil_Tf")
        else:
            append_log(System.subdir,"Finished: Tf_loop_iteration")
    System.error = error
    os.chdir(cwd)

def gmxcheck_subdirectories():
    """ Run gmxcheck on all traj.xtc files in subdirecories. """
    runs = glob("*/traj.xtc")
    dirs = [ x[:-9] for x in runs ]
    for subdir in dirs:
        run_gmxcheck(subdir)
    

def run_gmxcheck(dir):
    cwd = os.getcwd()
    os.chdir(dir)
    print "  Running gmxcheck on ",dir
    check = "gmxcheck -f traj.xtc"
    sb.call(check.split(),stdout=open("check.out","w"),stderr=open("check.err","w")) 
    
    error = "Fatal error"
    if (error in open("check.err","r").read()) or (error in open("check.out","r").read()):
        print "  FATAL ERROR in directory: ",dir
        print "  somethings wrong with Gromacs traj.xtc file. See %s/check.err" % dir
        #print open("check.err","r").read()
        error = 1
    else:
        error = 0
    os.chdir(cwd)
    return error

def determine_new_T_array():
    """ Find the temperatures which bracket the folding temperature.
        This takes the temperatures at which the average fraction of
        nonhelical contacts falls below 0.5 as bracketing the folding 
        temperature. A more complicated calculation is probably 
        needed for more complicated systems (proteins with intermediates)"""
    temps = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in temps ]
    temperatures.sort()
    lowerT = int(temperatures[0].split("_")[0])
    dT = int(temperatures[1].split("_")[0]) - lowerT
    upperT = int(temperatures[-1].split("_")[0])
    ## Temperatures bracket the folding temperature if the average 
    ## fraction of native contacts goes from greater than 0.5 to less
    ## than 0.5.
    for tdir in temperatures:
        Q = np.loadtxt(tdir+"/Qnhprob.dat")
        avgQ = np.mean(Q[len(Q)/2:])
        if avgQ > Q[0]/2:
            lowerT = int(tdir.split("_")[0])
        else:
            upperT = int(tdir.split("_")[0])
            break
        
    if dT == 2:
        ## Previous run was final iteration. Now WHAM needs to be
        ## done.
        newTi, newTf, newdT = 0,0,0
    else:
        ## Determine the finer grain T_array.
        newdT = int(float(dT)/5.)
        ## If newdT < 5 then just do finest temperature spacing. 
        if newdT < 5:
            newdT = 2
            midT = int(0.5*(float(lowerT)+upperT))
            newTi = midT - 20
            newTf = midT + 20
        else:
            newTi = lowerT + newdT
            newTf = upperT - newdT
    print "##DEBUGGING: New Ti, Tf, dT", newTi, newTf, newdT
    return newTi, newTf, newdT

def folding_temperature_loop(Model,System,append_log,new=False):
    """ The "folding temperature loop" is one of the several large-scale 
        logical structures in modelbuilder. It is entered anytime we want
        to determine the folding temperature. This could be when we have
        started a new project, refined the paramters, or returned to a 
        project in progress. The folding temperature loop successively 
        narrows in on the folding temperature."""

    cwd = os.getcwd()
    sub = System.path+"/"+System.subdir+"/"+System.Tf_active_directory
    #print sub  ## DEBUGGING
    if (not os.path.exists(sub)):
        os.mkdir(sub)
    ## Check to see if the folding temperature has been found. If yes, then continue.
    if (not os.path.exists(sub+"/Tf.txt")):
        os.chdir(sub)
        folding_temperature_loop_extension(Model,System,append_log,new=new)
    else:
        ## Folding temperature has been found. Continuing.
        pass
    os.chdir(cwd)

def folding_temperature_loop_extension(Model,System,append_log,new=False):
    """ This is for doing an additional loop in the Tf_loop. It either starts
        an initial temperature array or refines the temperature range according
        to previous data. """
    ## Check to see if a previous temperature range was used.
    if (not os.path.exists("T_array_last.txt")) or new:
        ## For initial exploration use very broad temperature increments.
        if System.initial_T_array != None:
            Ti = System.initial_T_array[0]
            Tf = System.initial_T_array[1]
            dT = System.initial_T_array[2]
        else:
            Ti = 50; Tf = 250; dT = 50
    else:
        ## Use previous range to determine new range. 
        Ti, Tf, dT = determine_new_T_array()
    print "  Running temperature array: T_initial = %.2f   T_final = %.2f   dT = %.2f " % (Ti,Tf,dT)
    run_temperature_array(Model,System,Ti,Tf,dT)
    append_log(System.subdir,"Submitting T_array iteration %d " % System.Tf_iteration)
    append_log(System.subdir,"  Ti = %d , Tf = %d , dT = %d" % (Ti, Tf, dT))
    append_log(System.subdir,"Starting: Tf_loop_iteration")

def start_next_Tf_loop_iteration(Model,System,append_log):
    """ To manually set the next temperature array."""

    Tf_choice = System.path+"/"+System.subdir+"/"+System.mutation_active_directory+"/Tf_choice.txt"
    Tf_guess = int(round(float(open(Tf_choice,"r").read()[:-1])))

    ## Update System counters
    System.Tf_iteration += 1
    System.Tf_active_directory = "Tf_"+str(System.Tf_iteration)
    System.mutation_iteration += 1
    System.mutation_active_directory = "Mut_"+str(System.mutation_iteration)

    cwd = os.getcwd()
    sub = System.path+"/"+ System.subdir+"/"+System.Tf_active_directory
    if os.path.exists(sub):
        print "ERROR!"
        print "  The next Tf iteration directory exists. "
        print "  exiting"
        raise SystemExit
    else:
        os.makedirs(sub)
    os.chdir(sub)
    Ti = int(Tf_guess - 20)
    Tf = int(Tf_guess + 20)
    dT = 2

    append_log(System.subdir,"Submitting T_array iteration %d" % System.Tf_iteration)
    append_log(System.subdir,"  Ti = %d , Tf = %d , dT = %d" % (Ti, Tf, dT))
    run_temperature_array(Model,System,Ti,Tf,dT)
    append_log(System.subdir,"Starting: Tf_loop_iteration")

    os.chdir(cwd)

def manually_add_temperature_array(Model,System,append_log,Ti,Tf,dT):
    """ To manually set the next temperature array."""
    cwd = os.getcwd()
    sub = System.path+"/"+ System.subdir+"/"+System.Tf_active_directory
    os.chdir(sub)
    append_log(System.subdir,"Submitting T_array iteration %d " % System.Tf_iteration)
    append_log(System.subdir,"  Ti = %d , Tf = %d , dT = %d" % (Ti, Tf, dT))
    run_temperature_array(Model,System,Ti,Tf,dT)
    append_log(System.subdir,"Starting: Tf_loop_iteration")

    os.chdir(cwd)

def run_equilibrium_simulations(Model,System,append_log):
    """ Run very long (equilibrium) simulations at the estimated folding 
        temperature."""

    if System.mutation_active_directory == '':
        System.mutation_active_directory = 'Mut_0'
        System.append_log("Creating mutational directory Mut_0" )
    
    cwd = os.getcwd()
    mutsub = System.path+"/"+System.subdir+"/"+System.mutation_active_directory
    Tfsub = System.path+"/"+System.subdir+"/"+System.Tf_active_directory
    Tf = open(Tfsub+"/Tf.txt","r").read().split()[0]

    System.append_log("Starting Equil_Tf")

    if not os.path.exists(mutsub):
        os.mkdir(mutsub)
    os.chdir(mutsub)
    T_string = ''
    for n in range(7):
        #T = "%.2f" % (float(Tf)+float(Tf)*(0.003*(n-1)))
        T = "%.2f" % (float(Tf)+float(Tf)*(0.003*(n-3)))
        for simnum in range(1,6):
            simpath = T+"_"+str(simnum)
            ## Only start the simulation if directory doesn't exist.
            if (not os.path.exists(simpath)):
                T_string += "%s\n" % simpath
                os.mkdir(simpath)
                os.chdir(simpath)
                System.append_log("  running T=%s" % simpath)
                #np.savetxt("Qref_cryst.dat",System.Qrefs[i],fmt="%1d",delimiter=" ")
                #print "Number of contacts: ", sum(sum(System.Qrefs[i]))
                print "    Running temperature ", T_string
                run_constant_temp(Model,System,float(T),nsteps=1000000000,walltime="60:00:00",queue="serial_long")
                os.chdir("..")
            else:
                ## Directory exists for this temperature: continue.
                continue

    open("T_array.txt","a").write(T_string)
    open("T_array_last.txt","w").write(T_string)
    append_log(System.subdir,"Starting: Equil_Tf")
    os.chdir(cwd)

def determine_walltime(Model):
    """ Estimate an efficient walltime."""
    Length = len(Model.Qref)
    ppn = "1"
    nsteps = "400000000"
    if Length < 60:
        walltime="12:00:00"
        queue="serial"
    else:
        if len(Model.Qref) > 160:
            nsteps = "600000000"
            if len(Model.Qref) > 250:
                nsteps = "800000000"
                walltime="72:00:00"
                ppn = "4"
            else:
                walltime="48:00:00"
                ppn = "2"
            queue="serial_long"
        else:
            walltime="24:00:00"
            queue="serial"
    return walltime, queue, ppn,nsteps

def run_temperature_array(Model,System,Ti,Tf,dT):
    """ Run many constant temperature runs over a range of temperatures to
        find the folding temperature. """

    System.append_log("Starting Tf_loop_iteration %d " % System.Tf_iteration)
    Temperatures = range(Ti,Tf+dT,dT)
    ## Run for longer if the protein is really big.
    walltime, queue, ppn, nsteps = determine_walltime(Model)

    T_string = ''
    for T in Temperatures:
        simpath = str(T)+"_0"
        ## Only start the simulation is directory doesn't exist.
        if (not os.path.exists(simpath)):
            T_string += "%d_0\n" % T
            os.mkdir(simpath)
            os.chdir(simpath)
            System.append_log("  running T=%d" % T)
            print "  Running temperature ", T
            run_constant_temp(Model,System,T,nsteps=nsteps,walltime=walltime,queue=queue,ppn=ppn)
            os.chdir("..")
        else:
            ## Directory exists for this temperature: continue.
            #print "##DEBUGGING: skipped ",T
            continue
    open("T_array.txt","a").write(T_string)
    open("T_array_last.txt","w").write(T_string)
    open("Ti_Tf_dT.txt","w").write("%d %d %d" % (Ti, Tf, dT))

def run_constant_temp(Model,System,T,nsteps="400000000",walltime="23:00:00",queue="serial",ppn="1"):
    ''' Start a constant temperature simulation with Gromacs. First it has
        to write the gromacs files stored in the System object, then it
        calls a function to submit the job.'''
    ## Loading and writing grompp.
    grompp_mdp = mdp.get_constant_temperature_mdp(Model,T,nsteps)
    open("grompp.mdp","w").write(grompp_mdp)

    ## Writing topology files.
    for filename in System.topology_files.iterkeys():
        #print "    Writing: ", filename    ## DEBUGGING
        open(filename,"w").write(System.topology_files[filename])

    ## Writing interaction tables. Writing native contact map.
    for m in range(len(Model.interaction_groups)):
        tablefile = "table_%s.xvg" % Model.interaction_groups[m]
        np.savetxt(tablefile,Model.tables[m],fmt="%16.15e",delimiter=" ")
    np.savetxt("table.xvg",Model.other_table,fmt="%16.15e",delimiter=" ")
    np.savetxt("Qref_cryst.dat",Model.Qref,fmt="%1d",delimiter=" ")
    ## Start simulation
    jobname = System.subdir+"_"+str(T)
    if Model.R_CD != None:
        jobname += "_Rcd_%.2f" % Model.R_CD
    if Model.dryrun == True:
        print "    Dryrun Success! Successfully saved simulation files." 
    else:
        submit_run(jobname,walltime=walltime,queue=queue,ppn=ppn)
    
def submit_run(jobname,walltime="23:00:00",queue="serial",ppn="1"):
    ''' Executes the constant temperature runs.'''

    prep_step1 = 'trjconv -f Native.pdb -o Native.xtc'
    prep_step2 = 'grompp -n index.ndx -c Native.pdb -maxwarn 1'
    prep_step3 = 'echo -e "System\\n" | trjconv -f Native.pdb -o conf.gro'
    #prep_step4 = 'grompp -n index.ndx -maxwarn 1'

    sb.call(prep_step1.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))
    sb.call(prep_step2.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))
    sb.call(prep_step3,shell=True,stdout=open("sim.out","w"),stderr=open("sim.err","w"))
    #sb.call(prep_step4.split())

    pbs_string = "#!/bin/bash \n"
    pbs_string +="#PBS -N %s \n" % jobname
    pbs_string +="#PBS -q %s \n" % queue
    pbs_string +="#PBS -l nodes=1:ppn=%s \n" % ppn
    pbs_string +="#PBS -l walltime=%s \n" % walltime
    pbs_string +="#PBS -V \n\n"
    pbs_string +="cd $PBS_O_WORKDIR\n"
    #pbs_string +="module purge\n"
    #pbs_string +="module load gromacs/4.6.5\n"
    pbs_string +="mdrun -s topol.tpr"

    open("run.pbs","w").write(pbs_string)
    qsub = "qsub run.pbs"
    sb.call(qsub.split(),stdout=open("sim.out","w"),stderr=open("sim.err","w"))

    rst_string = "#!/bin/bash \n"
    rst_string +="#PBS -N %s_rst \n" % jobname
    rst_string +="#PBS -q %s \n" % queue
    rst_string +="#PBS -l nodes=1:ppn=%s \n" % ppn
    rst_string +="#PBS -l walltime=%s \n" % walltime
    rst_string +="#PBS -V \n\n"
    rst_string +="cd $PBS_O_WORKDIR\n"
    #rst_string +="module purge\n"
    #rst_string +="module load gromacs/4.6.5\n"
    rst_string +="mdrun -s topol.tpr -cpi state.cpt"

    open("rst.pbs","w").write(rst_string)

if __name__ == '__main__':
    main()
