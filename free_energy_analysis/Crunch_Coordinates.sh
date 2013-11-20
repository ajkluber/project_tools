degen=$1
tempi=$2
tempstep=$3
tempf=$4
scale=$5

########################   Documentation   #######################
# Calculates the separate contributions from bonded and nonbonded terms.  It
# uses Gromacs g_energy and creates a data file where the contributions of the
# following terms to the potential energy. (The number for inputs to g_energy
# give different outputs depending on the Hamiltonian of the system.  So these
# input numbers are specific to the form of the Hamiltonitan used in DMC model
# 8-13-2013)
#
# output file is named energyterms.xvg, with rows corresponding to time and
# columns correspond to the different energetic terms, arranged as follows:
# 
# Columns
#  t 1 2 3 4 6
#
# Where
#  t = time in ps
#  1 = bonded energy
#  2 = angle energy
#  3 = dihedral energy
#  4 = nonboned energy
#  6 = total potential energy


WORKDIR=`pwd`
i=$1
for ((Temp=tempi;Temp<=tempf;Temp=Temp+tempstep));do
    #for ((i=0;i<degen;i++));do
    
        cd $WORKDIR
        cd ${Temp}_${i}

        # Use gromacs utility to calculate the rmsd with respect to native structure
        echo "Calculating rmsd..."
        echo -e "System\n" | trjconv -f traj.xtc -o whole.xtc -pbc whole
        echo -e "CA\nCA" | g_rms -f whole.xtc -s topol.tpr -o trajrmsd${Temp}_${i}.xvg -nomw -xvg none -n index.ndx
        cp trajrmsd${Temp}_${i}.xvg $WORKDIR/TrajectoryFiles
        cp trajrmsd${Temp}_${i}.xvg rmsd.xvg

        echo "Calculating radius of gyration..."
        echo "1" | g_gyrate -f traj.xtc -s topol.tpr -o radius_gyration.xvg -xvg none  #Radius of gyration
        awk '{print $1,$2}' radius_gyration.xvg > radius_cropped.xvg

        if [[ ! -e ener.edr ]]; then
            echo "ener.edr Does not exist for ${scale} T=${Temp}. Continuing..."
            continue
            exit 2
        else
            echo "Calculating energy terms..."
            cat > energterms.pbs << EOF
#!/bin/bash
#PBS -N Eterms_${scale}_${Temp}
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=00:02:00
#PBS -j oe
#PBS -V

cd \$PBS_O_WORKDIR
echo "1 2 3 4 6" | g_energy -f ener.edr -o energyterms -xvg none
EOF
            qsub energterms.pbs
        fi
        cd ..
    #done
done
echo "Finished."
