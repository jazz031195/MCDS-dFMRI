#!/bin/bash -l


process=""
loc=""
swell_perc=""
icvf=""
path=""

print_usage() {
  printf "Usage: 
  p : number of processes
  l : location (intra or extra)
  s : swell percentage
  i : icvf 
  P : path to folder "
}

while getopts n:p:l:s:i:P: opts; do
   case ${opts} in
      p) process=${OPTARG} ;;
      l) loc=${OPTARG} ;;
      s) swell_perc=${OPTARG} ;;
      i) icvf=${OPTARG} ;;
      P) path=${OPTARG} ;;
      *) print_usage
        exit 1 ;;
   esac
done

path_to_conf="$path/docs/conf_file_examples/gammaDistributedCylinders_run.conf"
path_to_model="$path/docs/conf_file_examples/gammaDistributedCylinders_mod.conf"

cp "$path_to_model" "$path_to_conf"
path_to_scheme="$path/docs/scheme_files/PGSE_sample_scheme_new.scheme"
path_to_axons="$path/instructions/demos/output/cylinders/Substrates/icvf_${icvf}${icvf}_swell_${swell_perc}_gamma_distributed_dyn_cylinder_list.txt"
path_to_exp="$path/instructions/demos/output/cylinders/icvf_${icvf}_swell_${swell_perc}_${loc}"

sed -i s:replace_exp_prefix:"${path_to_exp}":g "${path_to_conf}"
sed -i s:replace_scheme_file:"${path_to_scheme}":g "${path_to_conf}"
sed -i s:replace_axon_list:"${path_to_axons}":g "${path_to_conf}"
sed -i s:replace_process:"${process}":g "${path_to_conf}"
sed -i s:replace_loc:"${loc}":g "${path_to_conf}"