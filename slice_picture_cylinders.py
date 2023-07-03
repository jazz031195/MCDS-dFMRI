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
    n = 0
    with open(file ) as f:
        for e, line in enumerate(f.readlines()):
            if len(line.split(' ')) >4:
                n += 1
    return n

def get_first_cylinder_array(file, N):
    nbr_cylinders = N
    cylinder_array = np.zeros((nbr_cylinders,5))
    with open(file) as f:
        e = 0
        for line in f.readlines():

            if e< nbr_cylinders:
                
                if len(line.split(' ')) > 4:
                    
                    cylinder_array[e] = np.array([float(i) for i in line.split(' ')[:]])
                    e = e+1


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
        with open(cur_path + '/MCDC_Simulator_public-master/instructions/demos/output/cylinders/'+name+'.traj.txt') as f:
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
                z_min = 0.00
                z_max = 0.01

                indexes = [i for i,z in enumerate(zp) if z > z_min and z < z_max]

                xp = [xp[i] for i in indexes]
                yp = [yp[i] for i in indexes]
                zp = [zp[i] for i in indexes]

        return xp, yp, zp

def draw_cercles_with_trajectory(cylinder_array, size, file_name,swell = False):
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

def draw_cercles(cylinder_array, size, swell = False):
    r = []
    fig, ax = plt.subplots(ncols = 2) 

    print(cylinder_array)

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

def get_spheres_array(N,file):

    nbr_cylinders = N
    cylinder_array = np.zeros((nbr_cylinders,8))
    with open(file) as f:
        e = 0
        for line in f.readlines():

            if e< nbr_cylinders:
                
                if len(line.split(' ')) > 4:
                    
                    cylinder_array[e] = np.array([float(i) for i in line.split(' ')[:]])
                    e = e+1


    return cylinder_array

def main(name):

    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/Substrates/{name}.txt"
    size = get_size (file)
    axons = get_spheres_array(size,file)
    axons_slice = []

    for i,axon in enumerate(axons):
        print(axon)

        R = axon[-2]
        if axon[-1] == 1 :
            R = R*np.sqrt(1.01)


        axons_slice.append([axon[0], axon[1], R, axon[-1]])
    
    draw_cercles(np.array(axons_slice), size, swell = False)


def main2(name, traj_name):
    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/Substrates/{name}.txt"
    size = get_size (file)
    axons = get_spheres_array(size,file)
    axons_slice = []

    for i,axon in enumerate(axons):
        print(axon)

        R = axon[-2]
        if axon[-1] == 1 :
            R = R*np.sqrt(1.01)


        axons_slice.append([axon[0], axon[1], R, axon[-1]])
    
    draw_cercles_with_trajectory(np.array(axons_slice), size, traj_name, swell = False)

def check_distances(name):
    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/Substrates/{name}.txt"
    size = get_size (file)
    axons = get_spheres_array(size,file)

    for axon in enumerate(axons):
        axon = axon[1]
        x = axon[0]
        y = axon[1]
        r = axon[2]
        r = r*np.sqrt(1.01)
        for axon_ in enumerate(axons):
            axon_ = axon_[1]
            x_ = axon_[0]
            y_ = axon_[1]
            r_ = axon_[2]
            r_ = r_*np.sqrt(1.01)

            if ((x-x_)**2+(y-y_)**2< (r_+r)**2):
                return True
    return False
        

main2("icvf_0.7_swell_1.0_gamma_distributed_dyn_cylinder_list", "cyl_straight_try_extra")
#collides = check_distances("icvf_0.7_swell_1.0_gamma_distributed_dyn_cylinder_list")
#print(collides)