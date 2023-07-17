#!/bin/bash -l
# see make_conf_file.py to see the parameters to use

path="/scratch/ideriedm/Simulation/MCDS-dFMRI/MCDC_Simulator_public-master"

./MCDC_Simulator_public-master/compile.sh
./MCDC_Simulator_public-master/MC-DC_Simulator "${path}/docs/conf_file_examples/meshIntraInitialization_neuron.conf"
        