import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm
from mayavi import mlab

cur_path = os.getcwd()
conf_file_path = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"


def get_spheres_array(file):

    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 3):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                axons.append(spheres)
                spheres = []

    return axons

def draw_cercles(axons, swell = False):
    r = []

    [phi,theta] = np.mgrid[0:2*np.pi:12j,0:np.pi:12j]
    x_ = np.cos(phi)*np.sin(theta)
    y_ = np.sin(phi)*np.sin(theta)
    z_ = np.cos(theta)

    def plot_sphere(r,x,y,z):
        return mlab.mesh(r*x_+x, r*y_+y, r*z_+z)  
    for i in range(30):
        spheres = axons[i]
        for ii in range(len(spheres)):
            radius = spheres[ii][3]
            r.append(radius)
            x = spheres[ii][0]
            y = spheres[ii][1]
            z = spheres[ii][2]

        #circle1 = plt.Circle((x, y), radius, color='red')
        #ax[0].add_patch(circle1)

            plot_sphere(radius,x,y,z)

    mlab.show()


    if not swell:
        name = "slice.png"
    else:
        name = "slice_swell.png"

    plt.show()


def main():

    file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/_gamma_distributed_axon_list_.txt"

    axons = get_spheres_array(file)
    draw_cercles(axons)

main()


