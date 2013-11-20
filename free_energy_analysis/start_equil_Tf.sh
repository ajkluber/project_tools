WORKDIR=`pwd`
#degen=$1
#nsteps=400000000  ## 400000000 takes about 1 day.
nsteps=1000000000 ## 2.5x longer ~2.5days = 60hrs
degen=$1
scale=$2
jobprefix=$3

###############################      About      #############################
#
# For starting equilibrium simulations at the folding temperature, which is 
# assumed to be kept in file Tf_rounded.txt (Tf.txt). Starts a 
#
#
# Needs BeadBead.dat, Native.pdb, and Tf_rounded.txt. File Tf_rounded.txt
# should just be the folding temperature on the first line and nothing else.
###############################      About      #############################


if [[ ! -e Tf_rounded.txt ]]; then
    echo "Need Tf_rounded.txt. Exiting." 
    exit 2
else
    echo "Prepping..."
fi

#Tf=`cat Tf_rounded.txt`
Tfround=`cat Tf_rounded.txt`
Tf=`cat Tf.txt`

if [[ ! -d ${Tf}_${degen} ]]; then
    mkdir ${Tf}_${degen}
else
    echo "Directory exists. Exiting."
    exit 2
fi

if [[ -d Case_1 ]]; then
    if [[ -e Case_1/BeadBead.dat ]] && [[ -e Case_1/Native.pdb ]]; then
        cp Case_1/Native.pdb ${Tf}_${degen}/
        cp Case_1/BeadBead.dat ${Tf}_${degen}/
    fi
else
    echo "Native.pdb, BeadBead.dat, and t. Exiting"
fi

cd ${Tf}_${degen}
if [[ ! -e Native.xtc ]]; then
    echo "Starting new folding temperature run." >> Tfrun.log
    echo "Creating following files: index.ndx, Native.xtc, *.itp, topol.top." >> Tfrun.log
    echo -e "quit\n" | make_ndx -f Native.pdb
    trjconv -f Native.pdb -o Native.xtc 
    
    python $PROJECTS/dmc_model/Write_bond.py Native.pdb 1
    python $PROJECTS/dmc_model/Write_nonbond.py $scale
    echo "Copying files from ${Tf}_0: protein.itp, topol.top, grommp.mdp, table*.xvg" >> Tfrun.log
    cp ${PROJECTS}/dmc_model/gmx/grompp_constant_t.mdp ./grompp.mdp
    cp ${PROJECTS}/dmc_model/gmx/*.itp .
    cp ${PROJECTS}/dmc_model/gmx/topol.top .
    cp ${PROJECTS}/dmc_model/gmx/table*.xvg .
    grompp -n index.ndx -c Native.pdb -maxwarn 1
    echo -e "System\n" | trjconv -f Native.pdb -o conf.gro

    echo "Editing grommp.mdp to run ${nsteps} at T = ${Tf}" >> Tfrun.log
    sed -i 's/ref_t                    = 230/ref_t                    = '${Tf}'/g' grompp.mdp
    sed -i 's/gen_temp                 = 230/gen_temp                 = '${Tf}'/g' grompp.mdp
    sed -i 's/nsteps                   = 400000000/nsteps                   = '${nsteps}'/g' grompp.mdp
else
    echo "Restarting..."
fi

echo "Using grompp to generate topol.tpr" >> Tfrun.log
grompp -n index.ndx -maxwarn 1  &> initialize.log

echo "Running with pbs file: long_run.pbs " >> Tfrun.log
nthreads=1
jobname=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
cat > long_run.pbs <<EOF
#!/bin/bash
### Number of nodes and procs/node [DO NOT CHANGE]
#PBS -l nodes=1:ppn=${nthreads},walltime=60:00:00
###PBS -W group_list=pbc
#PBS -q serial_long
#PBS -V
### output files
#PBS -o out
#PBS -e err
### Job Name (max 15 chars.) [CHANGE]
#PBS -N ${jobprefix}_${jobname}Tf                                                      

cd \$PBS_O_WORKDIR
mdrun -nt 1                                                     
EOF

echo "Submitting with job ID:"
qsub long_run.pbs

