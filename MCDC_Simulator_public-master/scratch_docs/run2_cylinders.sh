#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# icvf
icvf="0.5";
concentration=1000000;

declare -a locs=( "intra" "extra"); 

declare -a swells=( "0.2"); 

path="/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
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
            if [ "$l" == "intra" ]; 
            then
                n=$(echo "scale=2; $concentration * $icvf" | bc)
		n=${n%.*}
            fi
            if [ "$l" == "extra" ]; 
            then
                new_icvf=$(echo "1 - $icvf" | bc)
                n=$(echo "scale=2; $concentration * $new_icvf" | bc)
		n=${n%.*}
            fi
            # run code to create axons list file
	    chmod u+x make_conf_file_cylinders.sh
            ./make_conf_file_cylinders.sh -i $icvf -l "$l" -p 48 -n "$n" -s "$s" -P "$path"
            ./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/gammaDistributedCylinders_run.conf"
        done
        
    done
done
  

