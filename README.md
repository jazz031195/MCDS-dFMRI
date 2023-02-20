# Introduction

This code allows for the creation of a phantom substrate and the simulation of random walkers (water molecules) in it. 

# Make configuration files

make_conf_file.py can be used to create a configuration file with desired parameters. 

- a : number of axons
- t : number of time points
- l : localisation (intra or extra)
- s : state (active or rest)
- i : icvf
- c : c2 for ODF (= 1 for no crossing fibers)
- x : creates a substrate if true or reads a substrate from a file if false

# Run simulation 
To run the simulation, modify the run.sh to adapt the parameters. Then simply ./run.sh.


 