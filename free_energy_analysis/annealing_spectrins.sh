#!/bin/bash

WORKDIR=`pwd`
setnumber=$1

###############################      About      #############################
#
# For starting annealing simulations for the mutated spectrins. This program
# needs to be started from the mutations directory
# 
###############################      About      #############################


if [[ ! -e NewBeadBead1.dat ]]; then
    echo "Need NewBeadBead1.dat. Exiting." 
    exit 2
elif [[ ! -e NewBeadBead2.dat ]]; then
    echo "Need NewBeadBead2.dat. Exiting."
    exit 2
elif [[ ! -e NewBeadBead3.dat ]]; then
    echo "Need NewBeadBead3.dat. Exiting."
    exit 2
else
    echo "Prepping..."
fi

if [[ $setnumber -eq 7 ]]; then
    sim_dir='seventh'
elif [[ $setnumber -eq 8 ]]; then
    sim_dir='eighth'
else
    echo "Please enter 7 or 8 (according to whether the data correspond to either directory)"
    exit 2
fi

for dir in r15 r16 r17; do 
    mkdir $dir
done

cp NewBeadBead1.dat ./r15/BeadBead.dat
cp NewBeadBead2.dat ./r16/BeadBead.dat
cp NewBeadBead3.dat ./r17/BeadBead.dat

for dir in r15 r16 r17; do
    
    cd $dir

    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.itp .
    rm nonbond_params.itp
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/index.ndx .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.pdb .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.xtc .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.gro .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.top .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.mdp .
    cp ~/projects/backup2/${sim_dir}-r15r16r17/${dir}-multi/*.xvg .
    python $PROJECTS/dmc_model/Write_nonbond.py 1
    echo "Using grompp to generate topol.tpr"
    grompp -n index.ndx -maxwarn 1 &> initialize.log
    cd ..
done

for dir in r15 r16 r16; do
    cd $dir
      
    echo "Running with pbs file: long_run.pbs on "`date` >> annealing_run.log
    nthreads=1
    jobname=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
    cat > run_${dir}.pbs <<EOF
#!/bin/bash
### Number of nodes and procs/node [DO NOT CHANGE]
#PBS -l nodes=1:ppn=${nthreads},walltime=24:00:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -o out
#PBS -e err
### Job Name (max 15 chars.) [CHANGE]
#PBS -N ${jobname}                                                      

cd \$PBS_O_WORKDIR
mdrun -nt 1                                                     
EOF

    echo "Submitting with job ID:"
    qsub run_${dir}.pbs
    cd ..
done