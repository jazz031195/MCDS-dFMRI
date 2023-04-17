#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

#c2
c2="1";
# nbr of axons
a="0";
# icvf
icvf="0.3";

#declare -a locs=("intra" "extra"); 
declare -a locs=("extra" "intra"); 

declare -a swells=("0.3" "0.4" "0.5" "0.6" "0.7"); 
#declare -a swells=("0.3"); 

cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/" || exit 
#./MCDC_Simulator_public-master/compile.sh
for s in "${swells[@]}";
do
    for l in "${locs[@]}";
    do
        echo "Starting creation of substrate ... "
        # run code to create axons list file
        python3 make_conf_file.py -a $a -i $icvf -c $c2 -x false -t 1000 -l $l -S $s -s active
        ./MCDC_Simulator_public-master/MC-DC_Simulator "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"
    done
    
done
  

