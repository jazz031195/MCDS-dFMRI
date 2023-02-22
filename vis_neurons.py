import matplotlib.pyplot as plt
import random
import numpy as np

with open('/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/_neurons_list.txt') as f:
    lines = f.readlines()
    fig   = plt.figure()
    ax    = plt.axes(projection='3d')

    for i in range(len(lines)):
        coords = lines[i].split(' ')
        # If it is the soma, plot it in any case
        if "Soma" in coords[0]:
            coords = lines[i-1].split(' ')
            # ax.scatter3D(float(coords[0]), float(coords[1]), float(coords[2]), s = 4/3*np.pi*float(coords[3])**3*1000)
            # draw sphere
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
            y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
            z = np.cos(v)*float(coords[3]) + float(coords[2])
            ax.plot_wireframe(x, y, z, color="r")
        if len(coords) > 2:
            # Plot only one sphere out of four for the dendrites (otherwise, too expensive)
            a = random.randint(1, 4)
            if a==1:
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                ax.plot_wireframe(x, y, z, color="b")
                ax.set_xlim([0, 0.1])
                ax.set_ylim([0, 0.1])
                ax.set_zlim([0, 0.1])
                # ax.scatter3D(float(coords[0]), float(coords[1]), float(coords[2]), s = 4/3*np.pi*float(coords[3])**3*1000)
    plt.show()