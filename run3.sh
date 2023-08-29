#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# icvf
vox="30";
icvf="0.50";
straight="false";
declare -a locs=("extra"); 


#declare -a time_steps=("1000" "5000" "10000" "30000"); 
declare -a time_steps=("5000" "10000" "30000"); 

declare -a conc=("10000000" "100000000" "1000000000"); 
#declare -a conc=("100000000"); 

#path="/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master"
./MCDC_Simulator_public-master/compile.sh

for l in "${locs[@]}";
    do
    for t in "${time_steps[@]}";
    do
        for c in "${conc[@]}";
        do
            for ((i=0;i<1;i++)); 
            do
                echo "Starting creation of substrate ... "
                echo " Location : $l"
                echo " Concentration : $c"
                echo " Time steps : $t"

                # run code to create axons list file
                chmod u+x make_conf_file.sh
                name="C_${c}_T_${t}_${l}"
                input="growth_icvf_0.50_cap_5_vox_30.swc"
                    
                ./make_conf_file.sh -i $icvf -l "$l" -p 15 -v "$vox" -s "0" -P "$path" -S "$straight" -c "$c" -t "$t" -n "$name" -I "$input"
                
                ./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/gammaDistributedAxons_run.conf"
                
            done
        done
        
    done
done
  

