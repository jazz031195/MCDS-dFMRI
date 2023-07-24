#!/bin/bash -l

# path="/scratch/ideriedm/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="MCDC_Simulator_public-master"

./MCDC_Simulator_public-master/compile.sh

declare -a N=( "100" "1000" "5000" "10000" "15000"); 

declare -a T=( "1000" "5000" "10000" "15000"); 

icvf="0.3";
p="48";
l="intra";

for n in "${N[@]}";
do
    for t in "${T[@]}";
    do
        for ((i=1;i<=5;i++)); 
        do
            echo "Starting creation of substrate ... "
            echo " N : $n"
            echo " T : $t"
            
            # run code to create neurons list file
	        chmod u+x make_conf_file.sh
            if [ "$n" -lt "$p" ];
            then
                ./make_conf_file.sh -i $icvf -l "$l" -p "$n" -n "$n" -t "$t" -s "$s" -P "$path"
            else
                ./make_conf_file.sh -i $icvf -l "$l" -p "$p" -n "$n" -t "$t" -s "$s" -P "$path"
            fi
            # ./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/Neurons_from_file_run.conf"
            ./build/MC-DC_Simulator "${path}/docs/conf_file_examples/Neurons_from_file_run.conf"
        done
        
    done
done

        