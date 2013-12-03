import subprocess as sb
import numpy as np
import os
import glob 



def analyze_temperature_array(System,i,append_log):
    ''' Analyze the '''
    cwd = os.getcwd()
    if System.error[i] == 0:
        sub = System.subdirs[i]+"/"+System.Tf_active_directory[i]
        os.chdir(cwd+"/"+sub)
        tempfile = open("T_array_last.txt","r").readlines()
        temperatures = [ temp[:-1] for temp in tempfile  ]

        ## TO DO: ADD CALCULATION OF NATIVE CONTACTS.
        if (not os.path.exists("Q.dat")):
            pass
            #calculate_Q()

        avgrmsd = np.zeros(len(temperatures),float)
        lowerT = 0
        upperT = 10000
        cwd2 = os.getcwd()
        for k in range(len(temperatures)):
            tdir = temperatures[k]
            #print cwd2+"/"+tdir ## DEBUGGING
            os.chdir(cwd2+"/"+tdir)
            if (not os.path.exists("rmsd.xvg")) or (not os.path.exists("radius_cropped.xvg")) \
               (not os.path.exists("energyterms.xvg")) or (not os.path.exists("phis.xvg")):
                crunch_coordinates()

            t,rmsd = np.loadtxt("rmsd.xvg",unpack=True)
            avgrmsd[k] = np.mean(rmsd[int(len(rmsd)/2):])
            os.chdir(cwd2)
        ## Find temperatures that bracket the folding temperature.
        print avgrmsd
        folded = list((avgrmsd < 0.4).astype(int))
        unfolded = list((avgrmsd > 0.8).astype(int))
        lowerT =  (temperatures[folded.index(0)])[:-2]
        upperT =  (temperatures[unfolded.index(1)])[:-2]
        print lowerT, upperT
        open("T_brackets.txt","w").write("%s %s" % (lowerT,upperT))
        os.chdir(cwd)
        append_log(System.subdirs[i],"Starting: Tf_loop_analysis")
    else:
        continue

def check_completion(System,i,append_log):

    pass

def crunch_coordinates():
    ''' Crunch the following reaction coordinates with Gromacs: rmsd, radius
        gyration, dihedrals, and potential energy terms.'''
    cmd1 = 'echo -e "CA\nCA" | g_rms -f traj.xtc -s topol.tpr -o rmsd.xvg -nomw -xvg none -n index.ndx'
    cmd2 = 'echo "1" | g_gyrate -f traj.xtc -s topol.tpr -o radius_gyration.xvg -xvg none'
    cmd3 = 'g_angle -n dihedrals.ndx -ov phis.xvg -all -type dihedral -xvg none'
    energy_pbs = "#!/bin/bash\n"
    energy_pbs +="#PBS -N Eterms_${scale}_${Temp}\n"
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

