#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

#c2
c2="1";
# nbr of axons
a="500";
# icvf
icvf="0.6";

cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/" || exit 
# compiles code
./MCDC_Simulator_public-master/compile.sh
for (( i=45; i<50; i++ ));
do
    echo "Starting creation of substrate ${i} ... "
    # run code to create axons list file
    python3 make_conf_file.py -a $a -i $icvf -c $c2 -x true -g $i
    ./MCDC_Simulator_public-master/src/main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

    rm "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/${a}_${i}_DWI.txt"

done
  

