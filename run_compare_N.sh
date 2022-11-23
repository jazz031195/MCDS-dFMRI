#!/bin/bash -l
declare -a locs=("extra");
declare -a states=("active");
#declare -a folders=("N_10_4" "N_5_10_4" "N_10_5" "N_5_10_5" "N_10_6");
#declare -a C=(10000 50000 100000 500000 1000000);
declare -a folders=("N_10_6");
declare -a N=(1000000);
#declare -a confs=("conf1" "conf2" "conf3");
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
  

