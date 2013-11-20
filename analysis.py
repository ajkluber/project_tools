import subprocess as sb
import numpy as np
import os
import glob 



def analyze_temperature_array(System,Tarrayfile="T_array.txt"):
    cwd = os.getcwd()
    for sub in System.subdirs:
        os.chdir(cwd+"/"+sub)
        tempfile = open(Tarrayfile,"r").readlines()
        temperatures = [ temp[:-1] for temp in tempfile  ]
        
        ## TO DO: ADD CALCULATION OF NATIVE CONTACTS.
        #calculate_Q()
        cwd2 = os.getcwd()
        for tdir in temperatures:
            print cwd2+"/"+tdir
            os.chdir(cwd2+"/"+tdir)
            try:
                check_files_exist(["rmsd.xvg","radius_cropped.xvg","energyterms.xvg","phis.xvg"])
            except:
                crunch_coordinates()
            os.chdir(cwd2)
        os.chdir(cwd)

def check_files_exist(files):
    ''' This will return an error if the files don't exist.'''
    for file in files:
        test = open(file,"r")
        test.close()

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

