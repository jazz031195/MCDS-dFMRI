#!/bin/bash -l

# makes configuration files which are the list of cylinders of the simulation

declare -a confs=("conf1" "conf2" "conf3");
loc="extra"
state="rest"
folder=""
simulator_path="build/MC-DC_Simulator"
mcdc_dir="MCDC_Simulator_public-master"
out_dir="$mcdc_dir/instructions/demos/output/dynamic_cylinders"
for conf in "${confs[@]}"
    do
    echo "configuration : $conf"    
    # modifies create_substrate.conf file and saves it to gammaDistributedCylinders.conf
    python3 make_conf_file.py -l "$loc" -n "10" -c "$conf" -s "$state" -f "$folder" -d "1"
    # runs simulator
    $simulator_path "$mcdc_dir/docs/conf_file_examples/gammaDistributedCylinders.conf"
    # removes useless DWI and information text
    rm "$out_dir/dyn_cylinder_DWI.txt"
    rm "$out_dir/dyn_cylinder_simulation_info.txt"
    # renames the cylinder list txt
    mv "$out_dir/dyn_cylinder_gamma_distributed_dyn_cylinder_list.txt" "$out_dir/${conf}.txt" 
    echo "configuration : $conf is done !"
    
done

