
# Script to start many scaling factors.
WORKDIR=`pwd`
tempi=200
tempf=500
tempstep=10
degen=0
testtemp=340
proteinname='Native'

module load openmpi/1.6.3-gcc

#for scale in 1 0.9 0.8 0.7 0.6 0.5; do  ## TESTING
for scale in 0; do  ## TESTING
#for scale in  1 0.9 0.8 0.7 0.6 0.5; do 
#for scale in  0.4 0.3 0.2 0.1; do 
    cd $WORKDIR
    if [[ -d $scale ]]; then
        echo "Analyzing $scale ..." 
    else
        echo "Directory ./$scale does not exist. Exiting"
        exit 2
    fi

    echo "Entering $scale/"
    cd $scale
    
    # Start by calculating the RMSD and potential energy. Should create
    # TrajectoryFiles/ directory.
    if [[ ! -e ${testtemp}_0/rmsd.xvg ]] || [[ ! -e ${testtemp}_0/energyterms.xvg ]] ; then
        echo "Performing rmsd, Rg, Energy_Terms calculations..."
        $PROJECTS/dmc_model/free_energy_analysis/Crunch_Coordinates.sh $degen $tempi $tempstep $tempf $scale &> crunch.log  
    else
        echo "RMSD, Rg, and potential energy already calculated. Continuing."
    fi

    #continue

    if [[ ! -e ${testtemp}_0/Qprob.dat ]] && [[ -e ${testtemp}_0/rmsd.xvg ]]; then
        echo "Calculating Native contacts reference matrix..."
        python $PROJECTS/dmc_model/Analysis/contact_calculator.py ref --method prob --refT 200

        echo "Calculating Q, Qh with 3 mpi jobs..."
        cat > qmpi1.pbs << EOF
#!/bin/bash
### Number of nodes and procs/node 
#PBS -l nodes=1:ppn=10,walltime=00:05:00
###PBS -W group_list=pbc
#PBS -q serial 
#PBS -V
### output files
#PBS -o Qout1
#PBS -e Qerr1
### Job Name (max 15 chars.) 
#PBS -N ${scale}_Qmpi.1

cd \$PBS_O_WORKDIR
mpirun -n 10 python /projects/cecilia/ajk8/dmc_model/Analysis/contact_calculator.py calc --Ti 200 --dT 10 --Tf 290 --method prob --refT 200
EOF
        job1=`qsub qmpi1.pbs`

        cat > qmpi2.pbs << EOF
#!/bin/bash
### Number of nodes and procs/node 
#PBS -l nodes=1:ppn=10,walltime=00:05:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -o Qout2
#PBS -e Qerr2
### Job Name (max 15 chars.) 
#PBS -N ${scale}_Qmpi.2

cd \$PBS_O_WORKDIR
mpirun -n 10 python /projects/cecilia/ajk8/dmc_model/Analysis/contact_calculator.py calc --Ti 300 --dT 10 --Tf 390 --method prob --refT 200
EOF

        job2=`qsub qmpi2.pbs`

        cat > qmpi3.pbs << EOF
#!/bin/bash
### Number of nodes and procs/node 
#PBS -l nodes=1:ppn=11,walltime=00:05:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -o Qout3
#PBS -e Qerr3
### Job Name (max 15 chars.) 
#PBS -N ${scale}_Qmpi.3

cd \$PBS_O_WORKDIR
mpirun -n 11 python /projects/cecilia/ajk8/dmc_model/Analysis/contact_calculator.py calc --Ti 400 --dT 10 --Tf 500 --method prob --refT 200
EOF

        job3=`qsub qmpi3.pbs`
        
        echo "Waiting for jobs to finish:"
        echo $job1 
        echo $job2 
        echo $job3 
    

        ####          WHAM  
        # Create data files for WHAM input.
        #if [[ ! -d wham ]]; then
        #    mkdir wham
        #fi
        #echo "Creating columns for WHAM..."
        #$PROJECTS/dmc_model/free_energy_analysis/Columns_for_Wham.sh $degen $tempi $tempstep $tempf gomodel R15A 
        #echo "Done"

        # Perform WHAM
        #Perform_WHAM.sh $degen $tempi $tempstep $tempf 
    fi


done
echo "Finished."

