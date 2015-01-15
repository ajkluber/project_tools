''' Functions for submitted analysis scripts as PBS jobs

Description:

    Utility to submit PBS jobs to calculate simulation observables
such as rmsd, Q, R_g, potential energy, dihedral angles, and number
of helical residues for Gromacs trajectories.
'''

import subprocess as sb
import numpy as np
import os 

def crunch_Q(name,contact_type,walltime="00:01:00",ppn="1",queue="serial"):
    ''' Submit PBS job to calculate sets of residue-residue contacts.

    Calculates contacts with structure-based models gromacs function g_kuh_sbm 
    '''

    contact_pbs = "#!/bin/bash\n"
    contact_pbs +="#PBS -N Q_"+name+"\n"
    contact_pbs +="#PBS -q %s\n" % queue
    contact_pbs +="#PBS -l nodes=1:ppn=%s\n" % ppn
    contact_pbs +="#PBS -l walltime=%s\n" % walltime
    contact_pbs +="#PBS -j oe\n"
    contact_pbs +="#PBS -V\n\n"
    contact_pbs +="cd $PBS_O_WORKDIR\n"
    if contact_type == "Gaussian":
        #contact_pbs +='g_kuh_sbm -s conf.gro -f traj.xtc -n contacts.ndx -noshortcut -abscut -cut 0.1 -qiformat list \n'
        contact_pbs +='g_kuh_sbm -s conf.gro -f traj.xtc -n native_contacts.ndx -o Q -noshortcut -abscut -cut 0.1\n'
    else:
        #contact_pbs +='g_kuh_sbm -s conf.gro -f traj.xtc -n contacts.ndx -noshortcut -noabscut -cut 0.2 -qiformat list \n'
        contact_pbs +='g_kuh_sbm -s conf.gro -f traj.xtc -n native_contacts.ndx -o Q -noshortcut -noabscut -cut 0.2\n'
    #contact_pbs +='mv qimap.out qimap.dat\n'
    contact_pbs +='mv Q.out Q.dat\n'
    open("contacts.pbs","w").write(contact_pbs)
    qsub = "qsub contacts.pbs"
    sb.call(qsub.split(),stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
    

def crunch_all(name,contact_type,walltime="00:01:00",ppn="1",n_tables=0):
    ''' Submit PBS job to calculate observables

    Calculates rmsd, radius gyration, dihedrals, and potential energy with 
    Gromacs utilities.
    '''
    analysis_pbs = '#!/bin/bash\n'
    analysis_pbs +='#PBS -N crnch_%s\n' % name
    analysis_pbs +='#PBS -q serial\n'
    analysis_pbs +='#PBS -l nodes=1:ppn=%s\n' % ppn
    analysis_pbs +='#PBS -l walltime=%s\n' % walltime
    analysis_pbs +='#PBS -j oe\n'
    analysis_pbs +='#PBS -V\n\n'
    analysis_pbs +='cd $PBS_O_WORKDIR\n'
    analysis_pbs +='g_energy -f ener.edr -o energyterms -xvg none << EOF\nBond\nAngle\nCoulomb-(SR)\nProper-Dih.\nPotential\nEOF'

    open("analysis.pbs","w").write(analysis_pbs)
    qsub = "qsub analysis.pbs"
    sb.call(qsub.split(),stdout=open("energyterms.out","w"),stderr=open("energyterms.err","w"))

def reorganize_qimap():
    ''' Parse a couple Q coordinates from qimap.dat '''
    G = np.loadtxt("radius_gyration.xvg")
    np.savetxt("Rg.xvg",G[:,0:2])

    contacts = np.loadtxt("contacts.dat")
    helical_contacts = (contacts[:,1] == contacts[:,0]+4)
    local_contacts = ((contacts[:,1] == contacts[:,0]+4).astype(int) + \
                  (contacts[:,1] == contacts[:,0]+5).astype(int) + \
                  (contacts[:,1] == contacts[:,0]+6).astype(int)).astype(bool)
    nonhelical_contacts = (helical_contacts == False)
    nonlocal_contacts = (local_contacts == False)

    qimap = np.loadtxt("qimap.dat")

    Qh = sum(qimap[:,helical_contacts].T)
    Qnh = sum(qimap[:,nonhelical_contacts].T)
    Qlocal = sum(qimap[:,local_contacts].T)
    Qnonlocal = sum(qimap[:,nonlocal_contacts].T)

    np.savetxt("Qh.dat",Qh)
    np.savetxt("Qnh.dat",Qnh)
    np.savetxt("Qlocal.dat",Qlocal)
    np.savetxt("Qnonlocal.dat",Qnonlocal)

def crunch_Nh(tol=40.):
    ''' Compute number of helical residues 

    A residue is in a helical conformation if it is making an i+4 contact and
    the two dihedral angles between it and its contact are within tol of 50 
    degrees (the 'helical' dihedral).
    '''
    Qh = np.loadtxt("Qhres.dat")
    phis = np.loadtxt("phis.xvg")
    phis = phis[:,2:]

    native_Nh = (phis[0,:] > (50-tol)).astype(int)*(phis[0,:] < (50+tol)).astype(int)
    num_nat_Nh = sum(native_Nh)

    dih1 = (phis[:,:-1] > (50-tol)).astype(int)*(phis[:,:-1] < (50+tol)).astype(int)
    dih2 = (phis[:,1:] > (50-tol)).astype(int)*(phis[:,1:] < (50+tol)).astype(int)
    Nh = dih1*dih2*Qh[:,:-4]

    #Nh_norm = sum(Nh.T)/float(num_nat_Nh)

    np.savetxt("Nh.dat",sum(Nh.T))
    np.savetxt("Nhres.dat",Nh)

