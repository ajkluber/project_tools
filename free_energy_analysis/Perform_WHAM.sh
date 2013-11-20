#!/bin/bash

WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4


if [[ -e input.wham ]]; then
    rm input.wham
fi

let number_of_steps=($last-$first)/$step+1
echo -e $number_of_steps>>input.wham 
for ((Temp=first;Temp<=last;Temp=Temp+step)); do
    for ((i=0;i<degen;i++)); do
        echo -e $Temp >> input.wham
        if [[ $degen==1 ]]
            then 
                echo -e "Input_for_WHAM_"$Temp"_0.dat" >> input.wham
            else
                echo -e "Input_for_WHAM_"$Temp"_"$i".dat" >>input.wham
        fi
    done
done

cat ranges_for_WHAM.wham >> input.wham

# A word should be said about the ranges_for_WHAM.wham file. The meaning of the
# rows is as follows:

# (Column number for first coordinate in Input_for_WHAM files) (second coordinate) (potential energy)
# (Lowest value) (step) (number of steps) ---> for first coordinate
# (Lowest value) (step) (number of steps) ---> for second coordinate
# (Lowest value) (step) (number of steps) ---> for potential energy
# Name of free energy file
# Name of heat capacity file
# (number of temperatures for which free energy & Cv are calculated) (first temperature) (step)
# Unknown parameter which is generally set to zero
# Ordinal number for the estimated Tf (meaning, which file number starting from the first one)
# Unknown parameter which is generally set to zero 

cp $PROJECTS/dmc_model/free_energy_analysis/WHAM .

cat > WHAM_calculations.pbs <<EOF
#!/bin/bash
#PBS -N WHAM
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=04:00:00
#PBS -j oe
#PBS -V                                                                         
  
echo "Job run on: "
cd \$PBS_O_WORKDIR

./WHAM<input.wham
EOF

qsub WHAM_calculations.pbs
