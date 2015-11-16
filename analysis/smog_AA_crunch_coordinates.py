""" Functions for submitted analysis scripts as PBS jobs

Description:

    Utility to submit PBS jobs to calculate simulation observables
such as rmsd, Q, R_g, potential energy, dihedral angles, and number
of helical residues for Gromacs trajectories.

"""

import subprocess as sb
import numpy as np
import os 

def crunch_Q(name,contact_type,walltime="00:03:00",ppn="1",queue="serial"):
    """ Submit SLURM job to calculate sets of residue-residue contacts.

    Calculates contacts with structure-based models gromacs function g_kuh_sbm 
    """

    contact_slurm = "#!/bin/bash\n"
    contact_slurm +="#SBATCH --job-name=Q_"+name+"\n"
    contact_slurm +="#SBATCH --partition=%s\n" % queue
    contact_slurm +="#SBATCH --nodes=1\n"
    contact_slurm +="#SBATCH --ntasks-per-node=%s\n" % ppn
    contact_slurm +="#SBATCH --time=%s\n" % walltime
    contact_slurm +="#SBATCH --export=ALL \n\n"
    contact_slurm +="cd $SLURM_SUBMIT_DIR\n"
    contact_slurm +="echo 'I ran on:'\n"
    contact_slurm +="cat $SLURM_JOB_NODELIST\n"
     
    contact_slurm +='g_kuh_sbm -s smog.gro -f traj.xtc -n smog_long.ndx -o Q -noshortcut -abscut -cut 0.1\n'
    
    contact_slurm +='mv Q.out Q.dat\n'
    open("contacts.slurm","w").write(contact_slurm)
    sbatch = "sbatch contacts.slurm"
    sb.call(sbatch.split(),stdout=open("contacts.out","w"),stderr=open("contacts.err","w"))
    

def crunch_all(name,contact_type,walltime="00:03:00",ppn="1",n_tables=0,queue="serial"):
    """ Submit SLURM job to calculate observables

    Calculates rmsd, radius gyration, dihedrals, and potential energy with 
    Gromacs utilities.
    """
    analysis_slurm = '#!/bin/bash\n'
    analysis_slurm +='#SBATCH --job-name=crnch_%s\n' % name
    analysis_slurm +='#SBATCH --partition=%s\n' % queue
    analysis_slurm +='#SBATCH --nodes=1\n'
    analysis_slurm +='#SBATCH --ntasks-per-node=%s\n' % ppn
    analysis_slurm +='#SBATCH --time=%s\n' % walltime
    analysis_slurm +='#SBATCH --export=ALL \n\n'
    analysis_slurm +='cd $SLURM_SUBMIT_DIR\n'
    analysis_slurm +='g_energy_sbm -f ener.edr -o energyterms -xvg none << EOF\nBond\nAngle\nCoulomb-(SR)\nProper-Dih.\nPotential\nEOF'

    open("analysis.slurm","w").write(analysis_slurm)
    sbatch = "sbatch analysis.slurm"
    sb.call(sbatch.split(),stdout=open("energyterms.out","w"),stderr=open("energyterms.err","w"))

def reorganize_qimap():
    """ Parse a couple Q coordinates from qimap.dat 

    DEPRECATED
    """
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
    """ Compute number of helical residues 

    A residue is in a helical conformation if it is making an i+4 contact and
    the two dihedral angles between it and its contact are within tol of 50 
    degrees (the 'helical' dihedral).

    DEPRECATED
    """
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

