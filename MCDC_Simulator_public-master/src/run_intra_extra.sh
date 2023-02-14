#!/bin/bash -l
cd "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/src" || exit 

./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_intra.conf"

./main "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons_extra.conf"
