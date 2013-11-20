
# Script to start many scaling factors.
WORKDIR=`pwd`
tempi=200
tempf=500
tempstep=10
degen=1
testtemp=340
proteinname='Native'

#for scale in 1; do  ## TESTING
for scale in 1 0.9 0.8 0.7 0.6 0.5; do  ## TESTING
#for scale in 0.4 0.3 0.2 0.1; do  ## TESTING
    cd $WORKDIR
    cd $scale

    if [[ ! -d wham ]]; then
        mkdir wham
    fi
    if [[ ! -e wham/Heat_rmsd_Rg.dat ]]; then
        echo "Prepping Heat Capcity run for scale $scale..."
        python $PROJECTS/dmc_model/Analysis/wham.py prep --output HeatCap --Ti $tempi --dT $tempstep --Tf $tempf --Tguess 370
        echo "Running Heat Capacity run for scale $scale..."
        python $PROJECTS/dmc_model/Analysis/wham.py run --output HeatCap --Ti $tempi --dT $tempstep --Tf $tempf --Tguess 370
    fi


    #continue

    #python $PROJECTS/dmc_model/Analysis/wham.py prep --output FreeEnergy --Ti $tempi --dT $tempstep --Tf $tempf
    if [[ -e wham/Heat_rmsd_Rg.dat ]]; then
        echo "Looks like Heat Capacity has been previously run..."
        ## 2D Free Energy prep
        if [[ ! -e wham/input_rmsd_Rg_F.wham ]] || [[ ! -e wham/input_Q_Rg_F.wham ]]; then
            echo "Prepping 2D Free Energy run for scale $scale..."
            python $PROJECTS/dmc_model/Analysis/wham.py prep --output FreeEnergy --Ti $tempi --dT $tempstep --Tf $tempf
        fi
        ## 1D Free Energy prep
        if [[ ! -e wham/input_rmsd_1D_F.wham ]] || [[ ! -e wham/input_Q_1D_F.wham ]]; then
            echo "Prepping 1D Free Energy run for scale $scale..."
            python $PROJECTS/dmc_model/Analysis/wham.py prep --output 1DFreeEnergy --Ti $tempi --dT $tempstep --Tf $tempf
        fi
    fi
    ## 2D Free Energy Plots
    if [[ -e wham/input_rmsd_Rg_F.wham ]] && [[ ! -e wham/FreeEnergy_Q_Rg.dat ]] && [[ -e wham/input_Qh_Qnh_F.wham ]]; then
        echo "Running 2D Free Energy run for scale $scale..."
        python $PROJECTS/dmc_model/Analysis/wham.py run --output FreeEnergy --Ti $tempi --dT $tempstep --Tf $tempf
    fi
    ## 1D Free Energy Plots
    if [[ -e wham/input_rmsd_1D_F.wham ]] && [[ ! -e wham/FreeEnergy_Q_1D.dat ]]; then
        echo "Running 1D Free Energy run for scale $scale..."
        python $PROJECTS/dmc_model/Analysis/wham.py run --output 1DFreeEnergy --Ti $tempi --dT $tempstep --Tf $tempf
    fi


done
echo "Finished."




