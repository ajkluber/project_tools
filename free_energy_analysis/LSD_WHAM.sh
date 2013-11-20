#!/bin/bash

WORKDIR=`pwd`

cat > LSD_WHAM.pbs <<EOF
#!/bin/bash
#PBS -N LSD_WHAM
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=04:00:00
#PBS -j oe
#PBS -V                                               
echo "Job run on: "
cd \$PBS_O_WORKDIR

./WHAM<inputLSD.wham
EOF

qsub LSD_WHAM.pbs
