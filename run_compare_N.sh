#!/bin/bash -l
declare -a locs=("intra" "extra");
declare -a states=("rest" "active");
declare -a folders=("N_10_4" "N_5_10_4" "N_10_5" "N_5_10_5" "N_10_6");
declare -a C=(10000 50000 100000 500000 1000000);
declare -a confs=("conf1" "conf2" "conf3");

for (( i=0; i<${#C[@]}; i++ ));
do
    for conf in "${confs[@]}"
    do
        for l in "${locs[@]}";
        do
            c=${C[i]};
            f=${folders[i]};
            for s in "${states[@]}"
            do
                echo "loc: $l, concentration : $c, conf : $conf, state : $s, folder : $f is starting"
                
                python3 make_conf_file.py -l "$l" -n "$c" -c "$conf" -s "$s" -f "$f"
                #path="/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/$s/$l/$f/$conf/conf.txt";
                cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/src" || exit 
                ./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"

                cd ..
                cd ..
                echo "loc: $l, concentration : $c, conf : $conf, state : $s, folder : $f is done"
                
            done
        done
    done
done
  

