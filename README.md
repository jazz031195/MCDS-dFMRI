# Introduction

This code allows for the creation of a phantom substrate and the simulation of random walkers (water molecules) in it. 

# Make configuration files

To create 3 configurations (conf1, conf2, conf3), run : ./create_configs.sh

In the MCDC_Simulator_public-master/instructions/demos/ouput folder, a dynamic_cylinder folder is created with the new conf1.txt, conf2.txt and conf3.txt files.

# Run simulation on a chosen configuration
To run the simulation, you must first make modifications in the run.sh file. run.sh has many loops so that the simulation runs for different locations (intra or extra), different states (active or rest), different folder names to save the output data (usually named after the number of molecules chosen), different number of walkers and different configurations (ex : con1, conf2, conf3).

Adjust run.sh to your needs and execute it with ./run.sh.

New folders will be created with your data organised in them.


 