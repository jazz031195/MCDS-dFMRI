#!/bin/bash -l

cd "/scratch/jnguyend/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master/src"
search_dir=''
src_files=''
for entry in *.cpp
do
  src_files="$src_files"" ""$entry"
done

#  sbatch "$entry"

command="g++ -O3 -std=c++11 -lpthread -std=c++0x -pthread -w -I.$src_files -o ../MC-DC_Simulator"


eval "$command"
