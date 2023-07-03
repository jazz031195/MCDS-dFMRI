import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


def get_axons_array(file):
    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 4):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                if (len(spheres)!= 0):
                    axons.append(spheres)
                spheres = []
            

    return axons

def change_icvf(file):
    axons = get_axons_array(file)
    with open(file) as f:
        for e,line in enumerate(f.readlines()):
            if e == 5:
                max_limit = float(line)

    AreaV = (max_limit) * (max_limit)*(max_limit)

    AreaC = 0

    for i, axon in enumerate(axons):

        ax_length = 0
        mean_rad= 0
        sum_rad = 0

        for j, sphere in enumerate(axon):
            if j >0:    
           
                diff = np.array([ sphere[0], sphere[1], sphere[2]] )-np.array([ axon[j-1][0], axon[j-1][1], axon[j-1][2]] )
                l = np.linalg.norm(diff)

                ax_length += l
                if sphere[4] == 1:
                    radius =  sphere[3]*np.sqrt(1.01)
                else:
                    radius = sphere[3]
                sum_rad += radius

        mean_rad = sum_rad/len(axon)


        AreaC += np.pi * mean_rad * mean_rad * ax_length


    print(AreaC / AreaV)

file = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/icvf_0.6_swell_0.1_gamma_distributed_axon_list.txt"
change_icvf(file)