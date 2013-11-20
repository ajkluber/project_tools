#!/bin/bash

WORKDIR=`pwd`
Option=$1
protein_name=$2


#protein=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
#protein_name=`echo ${protein} | sed s/"r"/"R"/g`


case $Option in
      
     'gomodel')

     cat > Sigma_ij_go_model.pbs <<EOF

#!/bin/bash
### Number of nodes and procs/node
#PBS -l nodes=1:ppn=1,walltime=00:30:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -j oe
### Job Name (max 15 chars.) [CHANGE]
#PBS -N Sigma_ij_go_model

echo "Job run on: "
cd \$PBS_O_WORKDIR
python $PROJECTS/dmc_model/Analysis/Go_model_contact_maker2.py ${protein_name}.pdb BeadBeadOLD.dat BeadBeadGO.dat
EOF

     qsub Sigma_ij_go_model.pbs
     ;;

    'native')
    
    cat > Sigma_ij_native.pbs <<EOF
#!/bin/bash
#PBS -N Sigma_ij_native
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=00:30:00
#PBS -j oe
#PBS -V  

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/GoModel_contact_maker3.py ${protein_name}.pdb BeadBeadOLD.dat BeadBeadnormal.dat
EOF

    qsub Sigma_ij_native.pbs
    ;;
    
    *)
    echo "Please select either gomodel or native options"
    ;;

esac
