#!/bin/bash -l

N=""
process=""
loc=""
swell_perc=""
icvf=""
path=""
straight=""

print_usage() {
  printf "Usage: 
  n : number of particles
  p : number of processes
  l : location (intra or extra)
  s : swell percentage
  i : icvf 
  P : path to folder 
  S : if axons are straight : true, else false "
}

while getopts n:p:l:s:i:P:S: opts; do
   case ${opts} in
      n) N=${OPTARG} ;;
      p) process=${OPTARG} ;;
      l) loc=${OPTARG} ;;
      s) swell_perc=${OPTARG} ;;
      i) icvf=${OPTARG} ;;
      P) path=${OPTARG} ;;
      S) straight=${OPTARG} ;;
      *) print_usage
        exit 1 ;;
   esac
done

path_to_conf="MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_run.conf"
path_to_model="MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_mod.conf"

cp "$path_to_model" "$path_to_conf"
path_to_scheme="$path/docs/scheme_files/PGSE_sample_scheme_new.scheme"

if [ "$straight" == "true" ]; then
  path_to_axons="$path/instructions/demos/output/axons/Substrates/icvf_${icvf}_swell_${swell_perc}_straight_gamma_distributed_axon_list.txt"
  path_to_exp="$path/instructions/demos/output/axons/icvf_${icvf}_swell_${swell_perc}_straight_${loc}"
  sed -i s:replace_tortuous:"false":g "${path_to_conf}"
else
  path_to_axons="$path/instructions/demos/output/axons/Substrates/icvf_${icvf}_swell_${swell_perc}_gamma_distributed_axon_list.txt"
  path_to_exp="$path/instructions/demos/output/axons/icvf_${icvf}_swell_${swell_perc}_${loc}"
  sed -i s:replace_tortuous:"true":g "${path_to_conf}"
fi

sed -i s:replace_exp_prefix:"${path_to_exp}":g "${path_to_conf}"
sed -i s:replace_N:"${N}":g "${path_to_conf}"
sed -i s:replace_scheme_file:"${path_to_scheme}":g "${path_to_conf}"
sed -i s:replace_axon_list:"${path_to_axons}":g "${path_to_conf}"
sed -i s:replace_process:"${process}":g "${path_to_conf}"
sed -i s:replace_loc:"${loc}":g "${path_to_conf}"