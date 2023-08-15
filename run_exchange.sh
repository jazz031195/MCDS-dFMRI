#!/bin/bash -l

# path="/scratch/ideriedm/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"
path="MCDC_Simulator_public-master"

declare -a N=("10000"); 

declare -a T=("10000"); 

icvf="0.3";
p="24";
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
            
            ./build/MC-DC_Simulator "${path}/docs/conf_file_examples/Neurons_from_file_intra.conf"
        done
        
    done
done