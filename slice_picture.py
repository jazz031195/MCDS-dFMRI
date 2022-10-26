import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm

cur_path = os.getcwd()
conf_file_path = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"
cyl_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dyn_cylinder_gamma_rest_extra_gamma_distributed_dyn_cylinder_list.txt"
cyl_swell_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dyn_cylinder_gamma_rest_extra_gamma_distributed_dyn_cylinder_swell_list.txt"
cyl_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/active/extra/N_10_5/dyn_cylinder_gamma_active_extra_gamma_distributed_dyn_cylinder_list.txt"

def get_T_N ():
    N = 0
    T = 0
    with open(conf_file_path) as f:
        for line in f.readlines():
            if line.split(' ')[0] == "N":
                N = int(line.split(' ')[1][:-1])
            elif line.split(' ')[0] == "T":
                T = int(line.split(' ')[1][:-1])
    return T, N

def get_nbr_cylinders ():
    nbr_cylinders = 0
    with open(conf_file_path) as f:
        for line in f.readlines():
            if line.split(' ')[0] == "num_dyn_cylinders":
                nbr_cylinders = int(line.split(' ')[1][:-1])
    return nbr_cylinders

def get_cylinder_array(swell = False):
    nbr_cylinders = 6000
    cylinder_array = np.zeros((nbr_cylinders,8))
    if swell :
        file = cyl_swell_file_path
    else :
        file = cyl_file_path
    with open(file) as f:
        e = 0
        for line in f.readlines():

            if e< nbr_cylinders:
                if len(line.split(' ')) > 5:
                    
                    cylinder_array[e] = np.array([float(i)*0.001 for i in line.split(' ')])
                    e = e+1
    print(e)
    return cylinder_array

def draw_cercles(cylinder_array, swell = False):
    fig, ax = plt.subplots() 
    for i in range(cylinder_array.shape[0]):
        radius = cylinder_array[i][-2]
        x = cylinder_array[i][0]
        y = cylinder_array[i][1]
        circle1 = plt.Circle((x, y), radius, color='r')
        ax.add_patch(circle1)
        
    plt.xlim([-0.005,0.12])
    plt.ylim([-0.005,0.12])
    plt.show()
    if not swell:
        name = "slice.png"
    else:
        name = "slice_swell.png"
    fig.savefig(name)

def main():

    swell = False
    cylinder_array = get_cylinder_array(swell)
    draw_cercles(cylinder_array, swell)

main()
