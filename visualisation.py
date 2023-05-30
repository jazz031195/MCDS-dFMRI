import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from scipy.linalg import norm
cur_path = os.getcwd()
fig, axs = plt.subplots(1, 3, figsize=(15, 5))
name = "test_"

def traj_data(name):
        with open(cur_path + '/MCDC_Simulator_public-master/instructions/demos/output/axons/'+name+'.traj.txt') as f:
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

def plot(xp, yp, zp):
    axs[0].plot(xp, yp, 'g.', markersize=1)
    axs[1].plot(xp, zp, 'g.', markersize=1)
    axs[2].plot(zp, yp, 'g.', markersize=1)
    axs[0].set_xlim([0, 0.1])
    axs[0].set_ylim([0, 0.1])
    axs[0].set_xlabel("X [mm]")
    axs[0].set_ylabel("Y [mm]")
    axs[1].set_xlim([0, 0.1])
    axs[1].set_ylim([0, 0.1])
    axs[1].set_xlabel("X [mm]")
    axs[1].set_ylabel("Z [mm]")
    axs[2].set_xlim([0, 0.1])
    axs[2].set_ylim([0, 0.1])
    axs[2].set_xlabel("Z [mm]")
    axs[2].set_ylabel("Y [mm]")

    plt.show()

xp, yp, zp = traj_data('_')
plot(xp, yp, zp)