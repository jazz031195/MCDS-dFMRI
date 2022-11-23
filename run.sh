#!/bin/bash -l

# location : "extra" or "intra"
declare -a locs=("extra" "intra"); 
#state : "active" or "rest"
declare -a states=("active" "rest");
# name of folder to save data to
declare -a folders=("N_10_6");
# number of molecules (later adjusted to location)
declare -a N=(1000000);
#name of configuration file
declare -a confs=("conf2");

icvf=0.7;

for (( i=0; i<${#N[@]}; i++ ));
do
    for conf in "${confs[@]}"
    do
        for l in "${locs[@]}";
        do

            n=${N[i]};
            if [[ $l == "extra" ]];then
                n=$(( (1-icvf)*n ));
            else
                n=$(( icvf*n ));
            fi
          
            f=${folders[i]};
            for s in "${states[@]}"
            do
                echo "loc: $l, num_walkers : $n, conf : $conf, state : $s, folder : $f is starting"
                
                python3 make_conf_file.py -l "$l" -n "$n" -c "$conf" -s "$s" -f "$f"
                #path="/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/$s/$l/$f/$conf/conf.txt";
                cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/src" || exit 
                ./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"

                cd ..
                cd ..
                echo "loc: $l, num_walkers : $n, conf : $conf, state : $s, folder : $f is done"
                
            done
        done
    done
done
  

