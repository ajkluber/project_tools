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

def run_temperature_array(Model,System,Ti=100,Tf=200,dT=10):
    ''' Run many constant temperature runs over a range of temperatures to
        find the folding temperature. '''

    Temperatures = range(Ti,Tf+dT,dT)
    for i in range(len(System.subdirs)):
        System.append_log("Starting T_array in directory: %s" % System.subdirs[i])
        T_string = ''
        for T in Temperatures:
            os.chdir(System.path)
            simpath = System.subdirs[i]+"/"+str(T)+"_0"
            ## Only start the simulation is directory doesn't exist.
            if (not os.path.exists(simpath)):
                T_string += "%d_0\n" % T
                os.mkdir(simpath)
                os.chdir(simpath)
                System.append_log("  running T=%d" % T)
                run_constant_temp(Model,System,T,i)
            else:
                ## Directory exists for this temperature: continue.
                continue
        open(System.path+"/"+System.subdirs[i]+"/T_array.txt","a").write(T_string)

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

    sb.call(prep_step1.split())
    sb.call(prep_step2.split())
    sb.call(prep_step3,shell=True)
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
    sb.call(qsub.split())
