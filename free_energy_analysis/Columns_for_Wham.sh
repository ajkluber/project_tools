#!/bin/bash

WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4
Option=$5
protein=$6

#protein=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
                                                                                                     
for ((Temp=first;Temp<=last;Temp=Temp+step));do
    for ((i=0;i<degen;i++));do

        if [[ ! -d $WORKDIR/wham ]]; then
            mkdir $WORKDIR/wham
        fi

        cd ${Temp}_${i}

	    case $Option in

		'dmc')
            awk '{print $2}' trajrmsd${Temp}_${i}.xvg > a.dat
            awk '{print $1}' Q_contacts_${Temp}.dat  > b.dat
            awk '{print $2}' radius_cropped.xvg > c.dat
            awk '{print $2}' pot_energy.xvg > d.dat
            awk '{print $1}' A_contacts_${Temp}.dat  > e.dat
            paste a.dat b.dat e.dat c.dat d.dat | awk '{print $1, $2, $3, $4, $5}'> Input_for_WHAM_${Temp}_${i}.dat
            rm e.dat
		    ;;

       'gomodel')
            awk '{print $2}' trajrmsd${Temp}_${i}.xvg > a.dat
            awk '{print $1}' Q_contacts_${Temp}.dat  > b.dat
            awk '{print $2}' radius_cropped.xvg > c.dat
            awk '{print $2}' pot_energy.xvg > d.dat
            paste a.dat b.dat c.dat d.dat | awk '{print $1, $2, $3, $4}'> Input_for_WHAM_${Temp}_${i}.dat
		    ;;

		'gohnh')
		    awk '{print $2}' trajrmsd${Temp}_${i}.xvg > a.dat
		    awk '{print $1}' Q_contacts_h_${Temp}.dat  > b.dat
		    awk '{print $2}' radius_cropped.xvg > c.dat
		    awk '{print $2}' pot_energy.xvg > d.dat
		    awk '{print $1}' Q_contacts_non_h_${Temp}.dat > e.dat
            paste a.dat b.dat e.dat c.dat d.dat | awk '{print $1, $2, $3, $4, $5}'> Input_for_WHAM_${Temp}_${i}.dat
            rm e.dat
		    ;;

		*)
		    echo "Please provide an option: dmc, gomodel, gohnh"
		    ;;

	    esac
	       

        tmp1=`echo $WORKDIR | sed s/enertris/'.'/`
        tmp2=`echo $tmp1 | cut -d . -f2`                                                                    
        cp Input_for_WHAM_${Temp}_${i}.dat ../wham/

        rm a.dat
        rm b.dat
        rm c.dat
        rm d.dat
	    cd $WORKDIR
    done
done
