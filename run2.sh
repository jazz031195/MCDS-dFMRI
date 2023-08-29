#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# icvf
vox="30";
icvf="0.50";
straight="false";
declare -a locs=( "intra" "extra"); 

declare -a swells=( "50" "70" "100"); 

#path="/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master"
./MCDC_Simulator_public-master/compile.sh

for ((i=1;i<=10;i++)); 
do
    for s in "${swells[@]}";
    do
        for l in "${locs[@]}";
        do
            echo "Starting creation of substrate ... "
            echo " Swelling percentage : $s"
            echo " Location : $l"

            # run code to create axons list file
	        chmod u+x make_conf_file.sh
            ./make_conf_file.sh -i $icvf -l "$l" -p 15 -v "$vox" -s "$s" -P "$path" -S "$straight"
            ./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/gammaDistributedAxons_run.conf"
        done
        
    done
done
  

