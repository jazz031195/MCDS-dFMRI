#!/bin/bash -l
# ./make_conf_file.sh -i 0.3 -l "intra" -p 48 -n 10 -t 100 -s 0 -P /Users/ideriedm/Documents/MCDS/MCDS-dFMRI/MCDC_Simulator_public-master 
N=""
T=""
process=""
loc=""
swell_perc=""
icvf=""
path=""

print_usage() {
  printf "Usage: 
  n : number of neurons
  t : number of timesteps
  p : number of processes
  l : location (intra or extra)
  s : swell percentage
  i : icvf 
  P : path to folder "
}

while getopts n:t:p:l:s:i:P: opts; do
   case ${opts} in
      n) N=${OPTARG} ;;
      t) T=${OPTARG} ;;
      p) process=${OPTARG} ;;
      l) loc=${OPTARG} ;;
      s) swell_perc=${OPTARG} ;;
      i) icvf=${OPTARG} ;;
      P) path=${OPTARG} ;;
      *) print_usage
        exit 1 ;;
   esac
done

path_to_conf="$path/docs/conf_file_examples/Neurons_from_file_run.conf"
path_to_model="$path/docs/conf_file_examples/Neurons_from_file_model.conf"

cp "$path_to_model" "$path_to_conf"
path_to_scheme="$path/docs/scheme_files/PGSE_sample_scheme_new.scheme"
path_to_neurons="$path/instructions/demos/output/neurons/intra/_neurons_list.txt"
path_to_exp="$path/instructions/demos/output/neurons/intra/n_${N}_T_${T}"

sed -i'' -e s:replace_exp_prefix:"${path_to_exp}":g "${path_to_conf}"
sed -i'' -e s:replace_N:"${N}":g "${path_to_conf}"
sed -i'' -e s:replace_T:"${T}":g "${path_to_conf}"
sed -i'' -e s:replace_scheme_file:"${path_to_scheme}":g "${path_to_conf}"
sed -i'' -e s:replace_neuron_list:"${path_to_neurons}":g "${path_to_conf}"
sed -i'' -e s:replace_process:"${process}":g "${path_to_conf}"
sed -i'' -e s:replace_loc:"${loc}":g "${path_to_conf}"
