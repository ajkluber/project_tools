import numpy as np
import subprocess as sb
import os

import mdp

'''
Created: November 17, 2013
Purpose:
    This module will be the library for submitting the simulation jobs.

Description:
    To be used as a library.

'''

def folding_temperature_loop(Model,System,append_log):
    ''' The "folding temperature loop" is one of the several large-scale 
        logical structures in modelbuilder. It is entered anytime we want
        to determine the folding temperature. This could be when we have
        started a new project, refined the paramters, or returned to a 
        project in progress. The folding temperature loop successively 
        narrows in on the folding temperature.'''

    for k in range(len(System.subdirs)):
        sub = System.path +"/"+ System.subdirs[k] +"/"+ System.Tf_active_directory[k]
        print sub  ## DEBUGGING
        if (not os.path.exists(sub)):
            os.mkdir(sub)
        ## Check to see if the folding temperature has been found. If yes, then continue.
        if (not os.path.exists(sub+"/Tf.txt")):
            ## Check to see if a previous temperature range was used.
            os.chdir(sub)
            folding_temperature_loop_extension(Model,System,k,append_log)
        else:
            ## Folding temperature has been found. Continuing.
            continue

def folding_temperature_loop_extension(Model,System,k,append_log):
    ''' This is for doing an additional loop in the Tf_loop. It either starts
        an initial temperature array or refines the temperature range according
        to previous data. '''
    if (not os.path.exists("Ti_Tf_dT.txt")):
        ## For initial exploration use very broad temperature increments.
        append_log(System.subdirs[k],"Submitting T_array iteration %d ; refinement %d" % \
                        (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
        append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (100, 300, 50))
        run_temperature_array(Model,System,k,Ti=100,Tf=300,dT=50)
        append_log(System.subdirs[k],"Starting: Tf_loop")
    else:
        ## Use previous range to determine new range. 
        ## DOESN'T WORK YET. USE ANALYSIS TO ESTIMATE THE BRACKETING TEMPS.
        Ti,Tf,dT = open("Ti_Tf_dT.txt","r").read().split()
        Ti = int(Ti); Tf = int(Tf); dT = int(dT)
        lowerT, upperT = open("T_brackets.txt","r").read().split()
        lowerT = int(lowerT); upperT = int(upperT)
        newdT = float(dT)/5.
        ## If new dT is less than 1 then don't do it.
        if newdT <= 3.:
            newdT = 1
            midT = int(0.5*(float(lowerT)+upperT))
            newTi = midT - 5
            newTf = midT + 5
        else:
            newTi = lowerT + newdT
            newTf = upperT - newdT

        run_temperature_array(Model,System,Ti=newTi,Tf=newTf,dT=newdT)
        append_log(System.subdirs[k],"Submitting T_array iteration %d ; refinement %d" % \
                        (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
        append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (newTi, newTf, dT))
        append_log(System.subdirs[k],"Starting: Tf_loop")
        System.Tf_refinements[k][System.Tf_iteration[k]] += 1

def run_temperature_array(Model,System,i,Ti=100,Tf=200,dT=10):
    ''' Run many constant temperature runs over a range of temperatures to
        find the folding temperature. '''

    Temperatures = range(Ti,Tf+dT,dT)
    System.append_log(System.subdirs[i],"Starting T_array in directory: %s" % System.subdirs[i])
    T_string = ''
    for T in Temperatures:
        simpath = str(T)+"_0"
        ## Only start the simulation is directory doesn't exist.
        if (not os.path.exists(simpath)):
            T_string += "%d_0\n" % T
            os.mkdir(simpath)
            os.chdir(simpath)
            System.append_log(System.subdirs[i],"  running T=%d" % T)
            run_constant_temp(Model,System,T,i)
            os.chdir("..")
        else:
            ## Directory exists for this temperature: continue.
            continue
    open("T_array.txt","a").write(T_string)
    open("Ti_Tf_dT.txt","w").write("%d %d %d" % (Ti, Tf, dT))

def run_constant_temp(Model,System,T,prot_num):
    ''' Start a constant temperature simulation with Gromacs. First it has
        to write the gromacs files stored in the System object, then it
        calls a function to submit the job.'''
    grompp_mdp = mdp.get_constant_temperature_mdp(Model,T)
    for filename in System.topology_files[prot_num].iterkeys():
        #print "Writing: ", filename    ## DEBUGGING
        open(filename,"w").write(System.topology_files[prot_num][filename])
    open("grompp.mdp","w").write(grompp_mdp)
    open("Native.pdb","w").write(System.native_pdbs[prot_num])
    for m in range(len(Model.interaction_groups)):
        tablefile = "table_%s.xvg" % Model.interaction_groups[m]
        np.savetxt(tablefile,Model.tables[m],fmt="%16.15e",delimiter=" ")
    np.savetxt("table.xvg",Model.other_table,fmt="%16.15e",delimiter=" ")
    ## Start simulation
    submit_run(System.subdirs[prot_num]+"_"+str(T))
    
def submit_run(jobname,walltime="23:00:00",queue="serial"):
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
    pbs_string +="### Number of nodes and procs/node \n"
    pbs_string +="#PBS -l nodes=1:ppn=1,walltime=%s \n" % walltime
    pbs_string +="###PBS -W group_list=pbc \n"
    pbs_string +="#PBS -q %s \n" % queue
    pbs_string +="#PBS -V \n"
    pbs_string +="### output files \n"
    pbs_string +="#PBS -o out \n"
    pbs_string +="#PBS -e err \n"
    pbs_string +="### Job Name (max 15 chars.) \n"
    pbs_string +="#PBS -N %s \n\n" % jobname
    pbs_string +="cd $PBS_O_WORKDIR\n"
    pbs_string +="mdrun -nt 1"

    open("run.pbs","w").write(pbs_string)
    qsub = "qsub run.pbs"
    sb.call(qsub.split(),stdout=open("sim.out","w"),stderr=open("sim.out"))
