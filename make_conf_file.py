# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys
import getopt

cur_path = os.getcwd()
conf_file = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders_old.conf"

def get_variables():
    C = None 
    conf = None 
    loc = None 
    state = None
    folder = None

    argv = sys.argv[1:]
  
    try:
        opts, args = getopt.getopt(argv, ":n:c:l:s:f:")
      
    except:
        print("Error")
  
    for opt, arg in opts:
        if opt in ['-n']:
            C = arg
        elif opt in ['-c']:
            conf = arg
        elif opt in ['-l']:
            loc = arg
        elif opt in ['-s']:
            state = arg
        elif opt in ['-f']:
            folder = arg
        
    return C, conf, loc, folder, state
        
      

def create_lines(C, conf, loc, folder, state):

    if state == "rest":
        active_state = "false"
    elif state == "active":
        active_state = "true"

    new_lines =[] 
    with open(conf_file) as f:
        for line in f.readlines():
            e = line.split(' ')
            if e[0] == "C":
                l = "C " + str(C)
            elif e[0] == "active_state":
                l = "active_state " + str(active_state)
            elif e[0] == "exp_prefix":
                l =  "exp_prefix /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder"
            elif e[0] == "dyn_cylinders_list":
                l = "dyn_cylinder_list /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(conf)+".txt"
            elif e[0] == "ini_walkers_pos":
                l = "ini_walkers_pos " + str(loc)
            else :
                l = line
            
            new_lines.append(l)
    print("created")
    return new_lines

def write_conf_file(lines, C, conf, loc, folder, state):

    name = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"

    f = open(name, "w")
    for l in lines :
        f.write(str(l))
        f.write("\n")
    f.close()

C, conf, loc, folder, state = get_variables()
lines = create_lines(C, conf, loc, folder, state)
write_conf_file(lines, C, conf, loc, folder, state)