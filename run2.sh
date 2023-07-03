#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# icvf
icvf_="0.7";
concentration=350000;
straight="true";
declare -a locs=( "intra" "extra"); 

declare -a swells=( "0.3" ); 

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
            if [ "$straight" == "true" ];
            then
            icvf="$(./get_icvf_axons.sh -f ${path}/instructions/demos/output/axons/Substrates/icvf_${icvf_}_swell_${s}_straight_gamma_distributed_axon_list.txt)"
            fi
            if [ "$straight" != "true" ];
            then
            icvf="$(./get_icvf_axons.sh -f ${path}/instructions/demos/output/axons/Substrates/icvf_${icvf_}_swell_${s}_gamma_distributed_axon_list.txt)"
            fi
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
	    chmod u+x make_conf_file.sh
            ./make_conf_file.sh -i $icvf_ -l "$l" -p 48 -n "$n" -s "$s" -P "$path" -S "$straight"
            ./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/gammaDistributedAxons_run.conf"
        done
        
    done
done
  

