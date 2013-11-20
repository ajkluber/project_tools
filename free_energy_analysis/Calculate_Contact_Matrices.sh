#!/bin/bash

WORKDIR=`pwd`

Temp=$1
Option=$2
protein_name=$3

#protein=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
#protein_name=`echo ${protein} | sed s/"r"/"R"/g`


cd ${Temp}_0

PBS_WORKDIR=`pwd`

case $Option in
      
'gomodel')

    
    cat > Q_matrix_go_model.pbs <<EOF
#!/bin/bash                                              
### Number of nodes and procs/node
#PBS -l nodes=1:ppn=1,walltime=6:00:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -o out
#PBS -e err
### Job Name (max 15 chars.) [CHANGE]
#PBS -N Q_matrix_go_model                                           

echo "Job run on: "
cd \$PBS_O_WORKDIR
python $PROJECTS/dmc_model/Analysis/Contact_maker_go_model.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix_${Temp}.dat           
EOF

    qsub Q_matrix_go_model.pbs
    ;;

'dmc')
    
    cat > Q_matrix.pbs <<EOF
#!/bin/bash
#PBS -N Q_matrix
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=00:30:00
#PBS -j oe
#PBS -V  

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Contact_maker.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix.dat A_matrix.dat trajrmsd${Temp}_0.xvg
EOF

    qsub Q_matrix.pbs
    ;;
    
'gohnh')
    
    cat > Q_matrix_go_h_non_h.pbs <<EOF
#!/bin/bash
#PBS -N Q_matrix_go_h_non_h.pbs
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=06:00:00
#PBS -j oe
#PBS -V

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Contact_maker_go_model_h_non_h.py ${protein_name}.pdb Q_matrix_h.dat Q_matrix_non_h.dat

EOF
    qsub Q_matrix_go_h_non_h.pbs
    ;;

'pdb')

    cat > Determine_Q_from_pdb.sh <<EOF
#!/bin/bash                 
#PBS -N Determine_Q_from_pdb
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=06:00:00
#PBS -j oe
#PBS -V                             

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Determine_Q_from_pdb_file.py ${protein_name}.pdb BeadBead.dat test.dat
EOF
    qsub Determine_Q_from_pdb.sh
	;;

*)
    echo "Please select either dmc, gomodel, gohnh or pdb options"
    ;;

esac



cd $WORKDIR
