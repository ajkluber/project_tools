import subprocess as sb
import numpy as np
import os 
import argparse

''' Library and command line utility for computing common reaction
    coordinates, such as rmsd, Q, R_g, potential energy, dihedral
    angles.'''


def crunch_Q(name,walltime="00:10:00"):
    ''' Calculate the fraction of native contacts, non-native contacts.'''

    contact_pbs = "#!/bin/bash\n"
    contact_pbs +="#PBS -N Q_"+name+"\n"
    contact_pbs +="#PBS -q serial\n"
    contact_pbs +="#PBS -l nodes=1:ppn=4,walltime=%s\n" % walltime
    contact_pbs +="#PBS -j oe\n"
    contact_pbs +="#PBS -V\n\n"
    contact_pbs +="cd $PBS_O_WORKDIR\n"
    contact_pbs +='python -m model_builder.analysis.contacts --calc  \n'
    open("contacts.pbs","w").write(contact_pbs)
    qsub = "qsub contacts.pbs"
    sb.call(qsub.split(),stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
    

def crunch_all(name,walltime="00:02:00"):
    ''' Crunch the following reaction coordinates with Gromacs: rmsd, radius
        gyration, dihedrals, and potential energy terms.'''
    cmd1 = 'echo -e "CA\nCA" | g_rms -f traj.xtc -s topol.tpr -o rmsd.xvg -nomw -xvg none -n index.ndx'
    cmd2 = 'echo "1" | g_gyrate -f traj.xtc -s topol.tpr -o radius_gyration.xvg -xvg none'
    cmd3 = 'g_angle -n dihedrals.ndx -ov phis.xvg -all -type dihedral -xvg none'
    energy_pbs = "#!/bin/bash\n"
    energy_pbs +="#PBS -N Eng_"+name+"\n"
    energy_pbs +="#PBS -q serial\n"
    energy_pbs +="#PBS -l nodes=1:ppn=4,walltime=%s\n" % walltime
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

def crunch_Nh(tol=40.):
    ''' Calculate number of helical residues per frame by the following criterion:
        if a residue is making an i+4 contact and the two dihedral angles between
        it and its contact'''
    Qh = np.loadtxt("Qhres.dat")
    phis = np.loadtxt("phis.xvg")
    phis = phis[:,2:]

    native_Nh = (phis[0,:] > (50-tol)).astype(int)*(phis[0,:] < (50+tol)).astype(int)
    num_nat_Nh = sum(native_Nh)

    dih1 = (phis[:,:-1] > (50-tol)).astype(int)*(phis[:,:-1] < (50+tol)).astype(int)
    dih2 = (phis[:,1:] > (50-tol)).astype(int)*(phis[:,1:] < (50+tol)).astype(int)
    Nh = dih1*dih2*Qh[:,4:]

    #Nh_norm = sum(Nh.T)/float(num_nat_Nh)

    np.savetxt("Nh.dat",sum(Nh.T))
    np.savetxt("Nhres.dat",Nh)


