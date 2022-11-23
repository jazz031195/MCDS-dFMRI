#!/bin/bash -l

# makes configuration files which are the list of cylinders of the simulation

declare -a confs=("conf1" "conf2" "conf3");
loc="extra"
state="rest"
folder=""
for conf in "${confs[@]}"
    do
    echo "configuration : $conf"    
    # modifies create_substrate.conf file and saves it to gammaDistributedCylinders.conf
    python3 make_conf_file.py -l "$loc" -n "10" -c "$conf" -s "$state" -f "$folder" -d "1"
    cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/src" || exit 
    # runs main
    ./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"
    # removes useless DWI and information text
    rm "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/dyn_cylinder_DWI.txt"
    rm "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/dyn_cylinder_simulation_info.txt"
    # renames the cylinder list txt
    mv "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/dyn_cylinder_gamma_distributed_dyn_cylinder_list.txt" "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/${conf}.txt" 
    cd ..
    cd ..
    echo "configuration : $conf is done !"
    
done

