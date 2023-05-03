import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm

cur_path = os.getcwd()
conf_file_path = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

def get_nbr_cylinders (file):
    nbr_cylinders = 0
    with open(file ) as f:
        for line in f.readlines():
            if "Axon" in line:
                nbr_cylinders += 1
    return nbr_cylinders

def get_size (file):
    with open(file ) as f:
        for e, line in enumerate(f.readlines()):
            if e == 5:
                return float(line)

def get_first_cylinder_array(file, N):
    nbr_cylinders = N
    cylinder_array = np.zeros((nbr_cylinders,5))
    prev_length = 0
    with open(file) as f:
        e = 0
        for line in f.readlines():

            if e< nbr_cylinders:
                
                if len(line.split(' ')) > 4 and prev_length < 4:
                    
                    
                    cylinder_array[e] = np.array([float(i) for i in line.split(' ')[:]])
                    e = e+1
            prev_length = len(line.split(' '))

    return cylinder_array


def get_last_cylinder_array(file, N):
    nbr_cylinders = N
    cylinder_array = np.zeros((nbr_cylinders,5))
    with open(file) as f:
        e = 0
        prev_line = ""
        line_nbr = 0
        for line in f.readlines():

            if e< nbr_cylinders:

                if (line_nbr>6 and len(line.split(' ')) <4 and len(prev_line.split(' '))>4):
                    
                    cylinder_array[e] = np.array([float(i) for i in prev_line.split(' ')[:]])
                    e = e+1
            prev_line= line
            line_nbr = line_nbr+1

    return cylinder_array

def traj_data(name):
        with open(cur_path + '/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/'+name+'.traj.txt') as f:
                xp = []
                yp = []
                zp = []
                i = 0
                for line in f.readlines():
                        if i%3 == 0:
                                xp.append(float(line))
                        elif i%3 == 1:
                                yp.append(float(line))
                        elif i%3 == 2:
                                zp.append(float(line))
                        i = i+1
                z_min = 0.001
                z_max = 0.079

                indexes = [i for i,z in enumerate(zp) if z > z_min and z < z_max]

                xp = [xp[i] for i in indexes]
                yp = [yp[i] for i in indexes]
                zp = [zp[i] for i in indexes]

        return xp, yp, zp

def draw_cercles(cylinder_array, size, file_name,swell = False):
    r = []
    fig, ax = plt.subplots(ncols = 2) 

    print(cylinder_array)

    xp, yp, zp = traj_data(file_name)

    for i in range(cylinder_array.shape[0]):
        radius = cylinder_array[i][-2]
        r.append(radius*1e3)
        x = cylinder_array[i][0]
        y = cylinder_array[i][1]
        swell = cylinder_array[i][-1]
        if (swell > 0):
            circle1 = plt.Circle((x, y), radius, color='red')
        else:
            circle1 = plt.Circle((x, y), radius, color='orange')
        ax[0].add_patch(circle1)
        
    ax[0].set_xlabel("[mm]")
    ax[0].set_ylabel("[mm]")
    ax[0].plot(xp, yp, 'g.', markersize=1)

    if not swell:
        name = "slice.png"
    else:
        name = "slice_swell.png"

    v = sns.histplot(data=pd.DataFrame({"radius[um]": r}),
            x="radius[um]", color='lightblue', ax = ax[1])
    ax[0].axis(xmin=0,xmax=size)
    ax[0].axis(ymin=0,ymax=size)
    plt.show()

    fig.savefig(name)

def get_spheres_array(file):

    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 3):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                if len(spheres) > 0:
                    axons.append(spheres)
                    spheres = []

    return axons

def main():

    file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/icvf_0.3_swell_0.7_gamma_distributed_axon_list.txt"
    N = get_nbr_cylinders (file )
    size = get_size (file)
    cylinder_array = get_first_cylinder_array(file, N)
    draw_cercles(cylinder_array, size)

    cylinder_array = get_last_cylinder_array(file, N)
    draw_cercles(cylinder_array, size)


def main2(name):
    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/{name}.txt"
    size = get_size (file)
    axons = get_spheres_array(file)
    axons_slice = []
    z = 0.01
    for axon in axons:

        index, value =  min(enumerate(list(np.array(axon).T[2])), key=lambda x: abs(x[1]-z))
        distance = np.abs(z-value)
        R = axon[index][3]
        if axon[index][4] == 1 :
            R = R*np.sqrt(1.01)
        if distance == 0:
            new_r = R
        else:
            if distance < R :
                new_r = np.sqrt(np.abs(R*R-distance*distance))
            else:
                continue
        axons_slice.append([axon[index][0], axon[index][1], z, new_r, axon[index][4]])
    
    draw_cercles(np.array(axons_slice), size, name, swell = False)
    
main2("icvf_0.3_swell_0.7_gamma_distributed_axon_list")