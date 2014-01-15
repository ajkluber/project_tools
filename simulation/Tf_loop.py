import numpy as np
import subprocess as sb
import os

import mdp

'''
Created: November 17, 2013
Purpose:
    This module will be the library for submitting the simulation jobs
for the Tf_loop (folding temperature loop).

Description:
    To be used as a library.

'''

def check_completion(System,i,append_log):
    ''' Checks to see if the previous Tf_loop simulation completed. First 
        checks the desired number of steps in the grompp.mdp file then 
        checks to see if md.log has recorded that number of steps.'''
    cwd = os.getcwd()
    sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
    os.chdir(cwd+"/"+sub)
    tempfile = open("T_array_last.txt","r").readlines()
    temperatures = [ temp[:-1] for temp in tempfile  ]
    print " Checking simulation completion..."
    for k in range(len(temperatures)):
        tdir = temperatures[k]
        error = 0
        ## Determine the number of steps for completed run.
        for line in open(tdir+"/"+"grompp.mdp","r"):
            if line[:6] == "nsteps":
                nsteps = int(line.split()[2]) + 1
                break    
        finish_line = "Statistics over " + str(nsteps)
        ## Check if md.log has finished the required number of steps.
        if finish_line in open(tdir+"/md.log","r").read():
            System.append_log(System.subdirs[i],"  %s finished." % tdir)
        else:
            print "    check %s. simulation did not finish."
            print "    Cannot continue with errors."
            System.append_log(System.subdirs[i],"  %s did not finish" % tdir)
            error = 1
    if error == 1:
        print "    Cannot continue with errors."
        pass 
    else:
        append_log(System.subdirs[i],"Finished: Tf_loop_iteration")
    System.error[i] = error
    os.chdir(cwd)

def determine_new_T_array():
    ''' Find the temperatures which bracket the folding temperature.
        This takes the temperatures at which the average fraction of
        native contacts falls below 0.5 as bracketing the folding 
        temperature. A more complicated calculation is probably 
        needed for more complicated systems (proteins with intermediates)'''
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
        Q = np.loadtxt(tdir+"/Qprob.dat")
        avgQ = np.mean(Q[len(Q)/2:])
        if avgQ > Q[0]/2:
            lowerT = int(tdir.split("_")[0])
        else:
            upperT = int(tdir.split("_")[0])
            break
        
    if dT == 1:
        ## Previous run was final iteration. Now WHAM needs to be
        ## done.
        newTi, newTf, newdT = 0,0,0
    else:
        ## Determine the finer grain T_array.
        newdT = int(float(dT)/5.)
        ## If newdT < 5 then just do finest temperature spacing. 
        if newdT < 5:
            newdT = 1
            midT = int(0.5*(float(lowerT)+upperT))
            newTi = midT - 5
            newTf = midT + 5
        else:
            newTi = lowerT + newdT
            newTf = upperT - newdT
    return newTi, newTf, newdT


def folding_temperature_loop(Model,System,k,append_log):
    ''' The "folding temperature loop" is one of the several large-scale 
        logical structures in modelbuilder. It is entered anytime we want
        to determine the folding temperature. This could be when we have
        started a new project, refined the paramters, or returned to a 
        project in progress. The folding temperature loop successively 
        narrows in on the folding temperature.'''

    cwd = os.getcwd()
    sub = System.path +"/"+ System.subdirs[k] +"/"+ System.Tf_active_directory[k]
    #print sub  ## DEBUGGING
    if (not os.path.exists(sub)):
        os.mkdir(sub)
    ## Check to see if the folding temperature has been found. If yes, then continue.
    if (not os.path.exists(sub+"/Tf.txt")):
        os.chdir(sub)
        folding_temperature_loop_extension(Model,System,k,append_log)
    else:
        ## Folding temperature has been found. Continuing.
        pass
    os.chdir(cwd)

def folding_temperature_loop_extension(Model,System,k,append_log):
    ''' This is for doing an additional loop in the Tf_loop. It either starts
        an initial temperature array or refines the temperature range according
        to previous data. '''
    ## Check to see if a previous temperature range was used.
    if (not os.path.exists("T_array_last.txt")):
        ## For initial exploration use very broad temperature increments.
        Ti = 50; Tf = 350; dT = 50
        append_log(System.subdirs[k],"Submitting T_array iteration %d ; refinement %d" % \
                        (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
        append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (Ti, Tf, dT))
        run_temperature_array(Model,System,k,Ti,Tf,dT)
        append_log(System.subdirs[k],"Starting: Tf_loop_iteration")
    else:
        ## Use previous range to determine new range. 
        newTi, newTf, newdT = determine_new_T_array()
        System.Tf_refinements[k][System.Tf_iteration[k]] += 1
        run_temperature_array(Model,System,k,newTi,newTf,newdT)
        append_log(System.subdirs[k],"Submitting T_array iteration %d ; refinement %d" % \
                        (System.Tf_iteration[k],System.Tf_refinements[k][System.Tf_iteration[k]]))
        append_log(System.subdirs[k],"  Ti = %d , Tf = %d , dT = %d" % (newTi, newTf, newdT))
        append_log(System.subdirs[k],"Starting: Tf_loop_iteration")

def run_temperature_array(Model,System,i,Ti,Tf,dT):
    ''' Run many constant temperature runs over a range of temperatures to
        find the folding temperature. '''

    Temperatures = range(Ti,Tf+dT,dT)
    System.append_log(System.subdirs[i],"Starting Tf_loop_iteration %d " % System.Tf_iteration[i])
    T_string = ''
    for T in Temperatures:
        simpath = str(T)+"_0"
        ## Only start the simulation is directory doesn't exist.
        if (not os.path.exists(simpath)):
            T_string += "%d_0\n" % T
            os.mkdir(simpath)
            os.chdir(simpath)
            System.append_log(System.subdirs[i],"  running T=%d" % T)
            run_constant_temp(Model,System,i,T)
            os.chdir("..")
        else:
            ## Directory exists for this temperature: continue.
            continue
    open("T_array.txt","a").write(T_string)
    open("T_array_last.txt","w").write(T_string)
    open("Ti_Tf_dT.txt","w").write("%d %d %d" % (Ti, Tf, dT))

def run_constant_temp(Model,System,k,T):
    ''' Start a constant temperature simulation with Gromacs. First it has
        to write the gromacs files stored in the System object, then it
        calls a function to submit the job.'''
    grompp_mdp = mdp.get_constant_temperature_mdp(Model,T)
    for filename in System.topology_files[k].iterkeys():
        #print "Writing: ", filename    ## DEBUGGING
        open(filename,"w").write(System.topology_files[k][filename])
    open("grompp.mdp","w").write(grompp_mdp)
    open("Native.pdb","w").write(System.native_pdbs[k])
    for m in range(len(Model.interaction_groups)):
        tablefile = "table_%s.xvg" % Model.interaction_groups[m]
        np.savetxt(tablefile,Model.tables[m],fmt="%16.15e",delimiter=" ")
    np.savetxt("table.xvg",Model.other_table,fmt="%16.15e",delimiter=" ")
    ## Start simulation
    submit_run(System.subdirs[k]+"_"+str(T))
    
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

def tail(f, window=1):
    ''' Quickly reads last line of long file.'''
    BUFSIZ = 1024
    f.seek(0, 2)
    bytes = f.tell()
    size = window
    block = -1
    data = []
    while size > 0 and bytes > 0:
        if (bytes - BUFSIZ > 0):
            # Seek back one whole BUFSIZ
            f.seek(block*BUFSIZ, 2)
            # read BUFFER
            data.append(f.read(BUFSIZ))
        else:
            # file too small, start from begining
            f.seek(0,0)
            # only read what was not read
            data.append(f.read(bytes))
        linesFound = data[-1].count('\n')
        size -= linesFound
        bytes -= BUFSIZ
        block -= 1
    return '\n'.join(''.join(data).splitlines()[-window:])

