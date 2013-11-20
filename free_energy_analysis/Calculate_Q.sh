#!/bin/bash

WORKDIR=`pwd`
degen=$1
first=$2
step=$3
last=$4
Option=$5
protein_name=$6

#protein=`pwd | sed s/"\/"/" "/g | awk '{print $NF}'`
#protein_name=`echo ${protein} | sed s/"r"/"R"/g`

for ((Temp=first;Temp<=last;Temp=Temp+step));do
    for ((i=0;i<degen;i++));do

        cd ${Temp}_${i}

        case $Option in

 		'dmc')
		
             cp ../Case_1/Q_matrix.dat .
             cp ../Case_1/A_matrix.dat .
             cat > Q_contacts.pbs <<EOF
#!/bin/bash
#PBS -N Q_dmc_${Temp}
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=02:00:00
#PBS -j oe
#PBS -V   

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Contact_reader.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix.dat A_matrix.dat Q_contacts_${Temp}.dat A_contacts_${Temp}.dat

EOF
            qsub Q_contacts.pbs
		;;
		
		'gomodel')
		 
             cp ../Case_1/Q_matrix.dat .
             cat > Q_contacts_gomodel.pbs <<EOF
#!/bin/bash
#PBS -N Q_Go_${Temp}
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=06:00:00
#PBS -j oe
#PBS -V   

echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Contact_reader_go_model.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix.dat Q_contacts_${Temp}.dat

EOF

             qsub Q_contacts_gomodel.pbs
	    
		 ;;

		 'gohnh')

              cp ../Case_1/Q_matrix_h.dat .
              cp ../Case_1/Q_matrix_non_h.dat .
              cat > Q_contacts_gohnh.pbs <<EOF
#!/bin/bash
#PBS -N Q_hnh_${Temp}
#PBS -q serial
#PBS -l nodes=1:ppn=1,walltime=06:00:00
#PBS -j oe
#PBS -V                                                                                        
echo "Job run on: "
cd \$PBS_O_WORKDIR

python $PROJECTS/dmc_model/Analysis/Contact_reader_go_model.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix_h.dat Q_contacts_h_${Temp}.dat
python $PROJECTS/dmc_model/Analysis/Contact_reader_go_model.py ${protein_name}.pdb BeadBead.dat whole.xtc Q_matrix_non_h.dat Q_contacts_non_h_${Temp}.dat

EOF
             qsub Q_contacts_gohnh.pbs

		 ;;

		 *)
		 
             echo "Please provide an option: dmc, gomodel, or gohnh"

		 ;;

	    esac
        cd $WORKDIR
    
    done
done

