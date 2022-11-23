# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
import sys
import getopt

cur_path = os.getcwd()

def get_variables():
    N = None 
    conf = None 
    loc = None 
    state = None
    folder = None
    create_substrate = 0

    argv = sys.argv[1:]
  
    try:
        opts, args = getopt.getopt(argv, ":n:c:l:s:f:d:")
      
    except:
        print("Error")
  
    for opt, arg in opts:
        if opt in ['-n']:
            N = arg
        elif opt in ['-c']:
            conf = arg
        elif opt in ['-l']:
            loc = arg
        elif opt in ['-s']:
            state = arg
        elif opt in ['-f']:
            folder = arg
        elif opt in ['-d']:
            create_substrate = int(arg)
        
        if create_substrate == 1:
            create_substrate = True
        else :
            create_substrate = False

        
        
    return N, conf, loc, folder, state, create_substrate
        
      

def create_lines(N, conf, loc, folder, state, create_substrate):
    if create_substrate:
        conf_file = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/create_substrate.conf"

    else:
        conf_file = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/read_from_substrate.conf"


    if state == "rest":
        active_state = "false"
    elif state == "active":
        active_state = "true"

    new_lines =[] 
    with open(conf_file) as f:
        for line in f.readlines():
            e = line.split(' ')
            if e[0] == "N":
                l = "N " + str(N)
            elif e[0] == "active_state":
                l = "active_state " + str(active_state)
            elif e[0] == "exp_prefix":
                if not create_substrate:
                    l = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)
                    if not os.path.exists(l):
                        os.makedirs(l)
                    l =  "exp_prefix " + l + "/dyn_cylinder"
                else :
                    l = "exp_prefix  /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/dyn_cylinder"
                    if not os.path.exists(l):
                        os.makedirs(l)
                    
            
            elif e[0] == "ini_walkers_pos":
                l = "ini_walkers_pos " + str(loc)
            
            elif e[0] == "dyn_cylinders_list":
                if not create_substrate:
                    l = "dyn_cylinder_list /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(conf)+".txt"
                else :
                    l = line
            else :
                l = line

            new_lines.append(l)

    return new_lines

def write_conf_file(lines):

    name = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"

    f = open(name, "w")
    for l in lines :
        f.write(str(l))
        f.write("\n")
    f.close()

N, conf, loc, folder, state, create_substrate = get_variables()
lines = create_lines(N, conf, loc, folder, state, create_substrate)
write_conf_file(lines)
