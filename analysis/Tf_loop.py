import subprocess as sb
import numpy as np
import os
import glob 

import contacts


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
                crunch_coordinates(System.subdirs[i]+"_"+tdir)
            os.chdir(cwd2)

        ## TO DO: ADD CALCULATION OF NATIVE CONTACTS.
        if (not os.path.exists(lowT+"/Qref_prob.dat")):
            os.chdir(lowT)
            Qref = contacts.probabilistic_reference()
            os.chdir(cwd2)
        else:
            Qref = np.loadtxt(lowT+"/Qref_prob.dat")

        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            np.savetxt("Qref_prob.dat",Qref,delimiter=" ",fmt="%1d")
            crunch_Q(System.subdirs[i]+"_"+tdir)

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


def crunch_Q(name):
    ''' Calculate the fraction of native contacts, non-native contacts.'''

    contact_pbs = "#!/bin/bash\n"
    contact_pbs +="#PBS -N Q_"+name+"\n"
    contact_pbs +="#PBS -q serial\n"
    contact_pbs +="#PBS -l nodes=1:ppn=1,walltime=00:05:00\n"
    contact_pbs +="#PBS -j oe\n"
    contact_pbs +="#PBS -V\n\n"
    contact_pbs +="cd $PBS_O_WORKDIR\n"
    contact_pbs +='python -m model_builder.analysis.contacts --calc  \n'
    open("contacts.pbs","w").write(contact_pbs)
    qsub = "qsub contacts.pbs"
    sb.call(qsub.split(),stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
    

def crunch_coordinates(name):
    ''' Crunch the following reaction coordinates with Gromacs: rmsd, radius
        gyration, dihedrals, and potential energy terms.'''
    cmd1 = 'echo -e "CA\nCA" | g_rms -f traj.xtc -s topol.tpr -o rmsd.xvg -nomw -xvg none -n index.ndx'
    cmd2 = 'echo "1" | g_gyrate -f traj.xtc -s topol.tpr -o radius_gyration.xvg -xvg none'
    cmd3 = 'g_angle -n dihedrals.ndx -ov phis.xvg -all -type dihedral -xvg none'
    energy_pbs = "#!/bin/bash\n"
    energy_pbs +="#PBS -N Eng_"+name+"\n"
    energy_pbs +="#PBS -q serial\n"
    energy_pbs +="#PBS -l nodes=1:ppn=1,walltime=00:02:00\n"
    energy_pbs +="#PBS -j oe\n"
    energy_pbs +="#PBS -V\n\n"
    energy_pbs +="cd $PBS_O_WORKDIR\n"
    energy_pbs +='echo "1 2 3 4 6" | g_energy -f ener.edr -o energyterms -xvg none\n'

    sb.call(cmd1,shell=True,stdout=open("rmsd.out","w"),stderr=open("rmsd.err","w"))
    sb.call(cmd2,shell=True,stdout=open("gyration.out","w"),stderr=open("gyration.err","w"))
    sb.call(cmd3,shell=True,stdout=open("dihedrals.out","w"),stderr=open("dihedrals.err","w"))
    G = np.loadtxt("radius_gyration.xvg")
    np.savetxt("radius_cropped.xvg",G[:,0:2])
    open("energyterms.pbs","w").write(energy_pbs)
    qsub = "qsub energyterms.pbs"
    sb.call(qsub.split(),stdout=open("energyterms.out","w"),stderr=open("energyterms.err","w"))

