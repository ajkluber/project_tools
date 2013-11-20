WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4

## NOW DEPRECATED
if [ ! -d $WORKDIR/TrajectoryFiles ]; then
    mkdir $WORKDIR/TrajectoryFiles
fi

for((Temp=first;Temp<=last;Temp=Temp+step));do
    for((i=0;i<degen;i++));do

        cd $WORKDIR
        cd ${Temp}_${i}

        # use gromacs utility to calculate the rmsd with respect to native structure
        echo -e "System\n" | trjconv -f traj.xtc -o whole.xtc -pbc whole
        echo -e "CA\nCA" | g_rms -f whole.xtc -s topol.tpr -o trajrmsd${Temp}_${i}.xvg -nomw -xvg none -n index.ndx
        cp trajrmsd${Temp}_${i}.xvg $WORKDIR/TrajectoryFiles
        cp trajrmsd${Temp}_${i}.xvg rmsd.xvg

        echo "1" | g_gyrate -f traj.xtc -s topol.tpr -o radius_gyration.xvg -xvg none  #Radius of gyration
        awk '{print $1,$2}' radius_gyration.xvg > radius_cropped.xvg
        #echo "6\n\n" | g_energy -f ener.edr -o pot_energy.xvg -xvg none
    
    done
done
