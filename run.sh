#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

# location : "extra" or "intra"
declare -a locs=("intra" "extra"); 
#state : "active" or "rest"
declare -a states=("rest" "active");
#number of time steps
T="1000";
#c2
c2="1";
# nbr of axons
a="50";
# icvf
declare -a icvf=("0.5");

cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/" || exit 
# compiles code
./MCDC_Simulator_public-master/compile.sh
for (( i=0; i<${#icvf[@]}; i++ ));
do
    echo "Starting creation of substrate with icvf of ${icvf[i]} ... "
    # run code to create axons list file
    python3 make_conf_file.py -a $a -i "${icvf[i]}" -c $c2 -x true
    ./MCDC_Simulator_public-master/src/main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

    rm "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/nbr_axons_{$a}_icvf_{${icvf[i]}}_DWI.txt"
    echo "Substrate with icvf of ${icvf[i]} is created !"
    for l in "${locs[@]}";
    do
        for s in "${states[@]}";
        do
            echo "Starting simulation in $l and in $s state ... "
            # run simulation on parameters of interest 
            python3 make_conf_file.py -a $a -i "${icvf[i]}" -c $c2 -s "$s" -l "$l" -t $T -x false
            ./MCDC_Simulator_public-master/src/main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"
            
            echo "Simulation in $l and in $s state is finished !"
        done
    done

done
  

