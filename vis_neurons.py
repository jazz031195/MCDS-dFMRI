import matplotlib.pyplot as plt
import random
import numpy as np
import os

wd = os.getcwd()

plot_3d = False
plot_traj = True

with open(wd + '/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/_neurons_list.txt') as f:
    lines = f.readlines()
    
    if plot_3d:
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim([0, 0.1])
        ax.set_ylim([0, 0.1])
        ax.set_zlim([0, 0.1])
        ax.set_xlabel('X [mm]')
        ax.set_ylabel('Y [mm]')
        ax.set_zlabel('Z [mm]')
    else:
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    for i in range(len(lines)):
        coords = lines[i].split(' ')
        # If it is the soma, plot it in any case
        if "Soma" in coords[0]:
            coords = lines[i-1].split(' ')
            # draw sphere
            u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
            x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
            y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
            z = np.cos(v)*float(coords[3]) + float(coords[2])
            if plot_3d:
                ax.plot_wireframe(x, y, z, color="r")
            else:
                axs[0].plot(x, y, color="r")
                axs[1].plot(x, z, color="r")
                axs[2].plot(z, y, color="r")
        elif len(coords) > 2:
            # Plot only one sphere out of four for the dendrites (otherwise, too expensive)
            a = random.randint(1, 10)
            if a==1:
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                if plot_3d:
                    ax.plot_wireframe(x, y, z, color="b")
                else:
                    axs[0].plot(x, y, color="b")
                    axs[1].plot(x, z, color="b")
                    axs[2].plot(z, y, color="b")
    if not plot_3d:
        distance_from_borders = 0.0007
        

    if plot_traj:
        with open(wd + '/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/.traj.txt') as f: 
            lines = f.readlines()
            xp = []
            yp = []
            zp = []
            for i in range(int(len(lines))):
                if i%3 == 0:
                    xp.append(float(lines[i]))
                elif i%3 == 1:
                    yp.append(float(lines[i]))
                elif i%3 == 2:
                    zp.append(float(lines[i]))
            if not plot_3d:
                axs[0].plot(xp, yp, 'g.', markersize=1)
                axs[1].plot(xp, zp, 'g.', markersize=1)
                axs[2].plot(zp, yp, 'g.', markersize=1)
                axs[0].set_xlim([0, 0.1])
                axs[0].set_ylim([0, 0.1])
                axs[0].set_xlabel("X [mm]")
                axs[0].set_ylabel("Y [mm]")
                axs[0].axvline(distance_from_borders, 0, 1)
                axs[0].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[0].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[0].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
                axs[1].set_xlim([0, 0.1])
                axs[1].set_ylim([0, 0.1])
                axs[1].set_xlabel("X [mm]")
                axs[1].set_ylabel("Z [mm]")
                axs[1].axvline(distance_from_borders, 0, 1)
                axs[1].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[1].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[1].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1)
                axs[2].set_xlim([0, 0.1])
                axs[2].set_ylim([0, 0.1])
                axs[2].set_xlabel("Z [mm]")
                axs[2].set_ylabel("Y [mm]")
                axs[2].axvline(distance_from_borders, 0, 1)
                axs[2].axvline(0.1-distance_from_borders, ymin=0, ymax=1)
                axs[2].axhline(distance_from_borders, xmin=0, xmax=1)
                axs[2].axhline(y=0.1-distance_from_borders, xmin=0, xmax=1) 
            else:
                ax.scatter(xp, yp, zp, color="g")

    plt.show()