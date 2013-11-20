WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4
scale=$5

###############################      Usage      #############################
#
# Many_TEMP <degen> <Ti> <dT> <Tf>
#
#
# Example: Start constant MD at temperature range 200-500 in steps of 10.
# Many_TEMP.sh 1 200 10 500
#
# Can start in empty directory with just a BeadBead.dat and Native.pdb file.
###############################      Usage      #############################


if [[ ! -e Native.pdb ]] && [[ ! -e BeadBead.dat ]]; then
    echo "Need Native.pdb and BeadBead.dat. Exiting."
    exit 2
else
    echo "Prepping Case_1 directory"
fi

if [[ ! -d Case_1 ]]; then
    mkdir Case_1
    cp Native.pdb Case_1/
    cp BeadBead.dat Case_1/
else
    if [[ ! -e Case_1/BeadBead.dat ]] || [[ ! -e Case_1/Native.pdb ]]; then
        cp Native.pdb Case_1/
        cp BeadBead.dat Case_1/
    fi
fi

cd Case_1
if [[ ! -e Native.xtc ]]; then
    echo "Creating following files: index.ndx, Native.xtc, *.itp, topol.top."
    echo -e "quit\n" | make_ndx -f Native.pdb
    trjconv -f Native.pdb -o Native.xtc 
    python $PROJECTS/dmc_model/Write_bond.py Native.pdb 1
    python $PROJECTS/dmc_model/Write_nonbond.py $scale
    cp $PROJECTS/dmc_model/gmx/protein.itp .
    cp $PROJECTS/dmc_model/gmx/topol.top .
    cp $PROJECTS/dmc_model/gmx/grompp_constant_t.mdp ./grompp.mdp
    cp $PROJECTS/dmc_model/gmx/table*.xvg .
    grompp -n index.ndx -c Native.pdb -maxwarn 1
    echo -e "System\n" | trjconv -f Native.pdb -o conf.gro
else
    echo "Case_1 prep directory exists. Copying..."
fi

cd $WORKDIR
# Prepare and Run Gromacs                                           
for ((Temp=first;Temp<=last;Temp=Temp+step));do
    for ((i=0;i<degen;i++));do

        if [[ ! -d ${Temp}_${i} ]] ; then
            #cp -r Case_1 ${Temp}_${i}
            mkdir ${Temp}_${i}
            cd ${Temp}_${i}
            echo "Editing grompp.mdp to correct temp: ${Temp}"
            sed -i 's/ref_t                    = 230/ref_t                    = '${Temp}'/g' grompp.mdp
            sed -i 's/gen_temp                 = 230/gen_temp                 = '${Temp}'/g' grompp.mdp
	    cd $WORKDIR
        fi
        cd ${Temp}_${i}

        echo "Using grompp on "
        grompp -n index.ndx -maxwarn 1  &> initialize.log
       
        nthreads=1
        jobname=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`

        cat > run.pbs <<EOF
#!/bin/bash
### Number of nodes and procs/node [DO NOT CHANGE]
#PBS -l nodes=1:ppn=${nthreads},walltime=23:00:00
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
        qsub run.pbs

        cd $WORKDIR
    done
done
