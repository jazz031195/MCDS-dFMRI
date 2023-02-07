import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm

cur_path = os.getcwd()
conf_file_path = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

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
            if line.split(' ')[0] == "num_axons":
                nbr_cylinders = int(line.split(' ')[1][:-1])
    return nbr_cylinders

def get_cylinder_array(file):
    nbr_cylinders = 50
    cylinder_array = np.zeros((nbr_cylinders,5))
    prev_length = 0
    with open(file) as f:
        e = 0
        for line in f.readlines():

            if e< nbr_cylinders:
                
                if len(line.split(' ')) > 4 and prev_length < 3:
                    
                    
                    cylinder_array[e] = np.array([float(i)*0.001 for i in line.split(' ')[:]])
                    e = e+1
            prev_length = len(line.split(' '))

    return cylinder_array

def draw_cercles(cylinder_array, swell = False):
    r = []
    fig, ax = plt.subplots(ncols = 2) 

    print(cylinder_array.shape)

    for i in range(cylinder_array.shape[0]):
        radius = cylinder_array[i][-2]
        r.append(radius*1000)
        x = cylinder_array[i][0]
        y = cylinder_array[i][1]
        swell = cylinder_array[i][-1]
        if (swell > 0):
            circle1 = plt.Circle((x, y), radius, color='red')
        else:
            circle1 = plt.Circle((x, y), radius, color='orange')
        ax[0].add_patch(circle1)
        
    #plt.xlabel("[mm]")
    #plt.ylabel("[mm]")

    if not swell:
        name = "slice.png"
    else:
        name = "slice_swell.png"

    v = sns.histplot(data=pd.DataFrame({"radius[um]": r}),
            x="radius[um]", kde=False, color='lightblue', ax = ax[1])

    plt.show()

    fig.savefig(name)

def main():

    file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/_gamma_distributed_axon_list.txt"

    cylinder_array = get_cylinder_array(file)
    draw_cercles(cylinder_array)

main()
