#!/bin/bash -l

# location : "extra" or "intra"
declare -a locs=("extra"); 
#state : "active" or "rest"
declare -a states=("rest");
# name of folder to save data to
declare -a folders=("N_10_6_");
# number of molecules (later adjusted to location)
declare -a N=(1000000);
#name of configuration file
declare -a confs=("conf1");

icvf=$(echo "0.7"|bc)

for (( i=0; i<${#N[@]}; i++ ));
do
    for conf in "${confs[@]}"
    do
        for l in "${locs[@]}";
        do

            n=${N[i]};
            if [[ $l == "extra" ]];then

                n=$(echo "(1-$icvf) * $n"|bc);

            else
                n=$(echo "$icvf * $n"|bc);

            fi
            n=${n%.*};

          
            f=${folders[i]};
            for s in "${states[@]}"
            do
                echo "loc: $l, num_walkers : $n, conf : $conf, state : $s, folder : $f is starting"
                
                python3 make_conf_file.py -l "$l" -n "$n" -c "$conf" -s "$s" -f "$f"
                cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/src" || exit 
                ./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"

                cd ..
                cd ..
                echo "loc: $l, num_walkers : $n, conf : $conf, state : $s, folder : $f is done"
                
            done
        done
    done
done
  

