WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4

###############################      Usage      #############################
# Used to restart runs that didn't make it to the number of steps specified in
# the topol.tpr file. The number of steps a simulation finished can be found in
# the md.log file (to get the number of completed frames use: grep "Statistics
# over" md.log). 
# 
#
#
# restart_runs.sh <degen> <Ti> <dT> <Tf>
#
#
# restart_runs.sh 1 200 10 500
#
# Restarts runs using the Gromacs saved checkpoint file. Info found at:
# http://www.gromacs.org/Documentation/How-tos/Doing_Restarts
###############################      Usage      #############################

# restart Gromacs                                           
for ((Temp=first;Temp<=last;Temp=Temp+step));do
    for ((i=0;i<degen;i++));do


	    cd $WORKDIR
        cd ${Temp}_${i}

        if [[ ! -e topol.tpr ]] || [[ ! -e state.cpt ]]; then
            echo "Topology file topol.tpr or checkpoint file state.cpt do not exist! Exiting."
            exit 2
        else
            echo "Restarting job with saved checkpoint file state.cpt"
        fi

        nthreads=1
        jobname=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`

        echo "Submitting job: "
        cat > restart_run.pbs <<EOF
#!/bin/bash
### Number of nodes and procs/node [DO NOT CHANGE]
#PBS -l nodes=1:ppn=${nthreads},walltime=23:00:00
###PBS -W group_list=pbc
#PBS -q serial
#PBS -V
### output files
#PBS -o rstout
#PBS -e rsterr
### Job Name (max 15 chars.) [CHANGE]
#PBS -N rst${jobname}                                                              

cd \$PBS_O_WORKDIR
mdrun -nt 1 -s topol.tpr -cpi state.cpt                                     
EOF

        qsub restart_run.pbs

    done
done
