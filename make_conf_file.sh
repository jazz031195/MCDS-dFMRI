#!/bin/bash -l

vox=""
process=""
loc=""
swell_perc=""
icvf=""
path=""
straight=""
concentration=""
time_steps=""
name=""
input=""

print_usage() {
  printf "Usage: 
  t : number of time steps
  c : concentration per mm3
  p : number of processes
  l : location (intra or extra)
  s : swell percentage
  i : icvf 
  P : path to folder 
  S : if axons are straight : true, else false 
  v : voxel size 
  n : name of file 
  I : name of input swc file "
}

while getopts t:c:n:p:l:s:i:P:S:v:I: opts; do
   case ${opts} in
      t) time_steps=${OPTARG} ;;
      c) concentration=${OPTARG} ;;
      p) process=${OPTARG} ;;
      l) loc=${OPTARG} ;;
      s) swell_perc=${OPTARG} ;;
      i) icvf=${OPTARG} ;;
      P) path=${OPTARG} ;;
      S) straight=${OPTARG} ;;
      v) vox=${OPTARG} ;;
      n) name=${OPTARG} ;;
      I) input=${OPTARG} ;;
      *) print_usage
        exit 1 ;;
   esac
done

path_to_conf="MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_run.conf"
path_to_model="MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_mod.conf"

cp "$path_to_model" "$path_to_conf"
path_to_scheme="$path/docs/scheme_files/PGSE_sample_scheme_new.scheme"

if [ "$straight" == "true" ]; then
  #path_to_axons="$path/instructions/demos/output/axons/straight/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}.swc"
  #path_to_exp="$path/instructions/demos/output/axons/straight/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}_${loc}"
  sed -i s:replace_tortuous:"false":g "${path_to_conf}"
else
  #path_to_axons="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}.swc"
  #path_to_exp="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/icvf_${icvf}_vox_${vox}_swell_${swell_perc}_${loc}"
  sed -i s:replace_tortuous:"true":g "${path_to_conf}"
fi
path_to_exp="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/${name}"
path_to_axons="$path/instructions/demos/output/axons/tortuous/icvf_${icvf}/${input}"
sed -i s:replace_exp_prefix:"${path_to_exp}":g "${path_to_conf}"
sed -i s:replace_scheme_file:"${path_to_scheme}":g "${path_to_conf}"
sed -i s:replace_axon_list:"${path_to_axons}":g "${path_to_conf}"
sed -i s:replace_process:"${process}":g "${path_to_conf}"
sed -i s:replace_loc:"${loc}":g "${path_to_conf}"
sed -i s:replace_concentration:"${concentration}":g "${path_to_conf}"
sed -i s:replace_timesteps:"${time_steps}":g "${path_to_conf}"