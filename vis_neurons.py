import matplotlib.pyplot as plt
import random
import numpy as np
import os
import plotly.graph_objs as go
wd = os.getcwd()
plot_3d = False
plot_traj = True
projection = True
z_slice = [0.02, 0.04, 0.06, 0.08]
max_lim = 0
neuron_file = wd + '/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/_rep_00_neurons_list.txt'
traj_file = wd + '/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/_rep_00.traj.txt'
with open(neuron_file) as f:
    lines = f.readlines()
    if plot_3d:
        plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlabel('X [mm]')
        ax.set_ylabel('Y [mm]')
        ax.set_zlabel('Z [mm]')
    elif projection:
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    else:
        fig, axs = plt.subplots(1, len(z_slice), figsize=(15, 5))
        for k in range(len(z_slice)):
            axs[k].set_xlabel("X [mm]")
            axs[k].set_ylabel("Y [mm]")
    # lines_plot = []
    for i in range(len(lines)):
        for j, slice_ in enumerate(z_slice):
            coords = lines[i].split(' ')
            # If it is the soma, plot it in any case
            if "Soma" in coords[0]:
                coords = lines[i-1].split(' ')
                coords = [float(coord) for coord in coords]
                if max_lim < coords[0]:
                    max_lim = coords[0]
                # draw sphere
                u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                z = np.cos(v)*float(coords[3]) + float(coords[2])
                if plot_3d:
                    ax.plot_wireframe(x, y, z, color="r")
                    # # Creating the plot
                    # line_marker = dict(color='red', width=2)
                    # for l, m, n in zip(x, y, z):
                    #     lines_plot.append(go.Scatter3d(x=l, y=m, z=n, mode='lines', line=line_marker))
                elif projection:
                    axs[0].plot(x, y, color="r")
                    axs[1].plot(x, z, color="r")
                    axs[2].plot(z, y, color="r")
                else:
                    if ((coords[2] - coords[3]) < slice_) and ((coords[2] + coords[3]) > slice_):
                        idx = ((x - coords[0])**2 + (y - coords[1])**2 < coords[3]**2 - (slice_ - coords[2])**2)
                        x = x[idx]
                        y = y[idx]
                        x = np.squeeze(x.reshape(1, -1))
                        y = np.squeeze(y.reshape(1, -1))
                        z = np.ones_like(x)*slice_
                        if z[0] == slice_:
                            axs[j].plot(x, y, '.',color="darkgray")
                    else:
                        z = np.squeeze(z.reshape(1, -1))
            elif len(coords) > 2:
                coords = [float(coord) for coord in coords]
                if max_lim < coords[0]:
                    max_lim = coords[0]
                # Plot only one sphere out of four for the dendrites (otherwise, too expensive)
                a = random.randint(1, 50)
                if a==1:
                    # draw sphere
                    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
                    x = np.cos(u)*np.sin(v)*float(coords[3]) + float(coords[0])
                    y = np.sin(u)*np.sin(v)*float(coords[3]) + float(coords[1])
                    z = np.cos(v)*float(coords[3]) + float(coords[2])
                    if plot_3d:
                        ax.plot_wireframe(x, y, z, color="r")
                        # # Creating the plot
                        # line_marker = dict(color='blue', width=2)
                        # for l, m, n in zip(x, y, z):
                        #     lines_plot.append(go.Scatter3d(x=l, y=m, z=n, mode='lines', line=line_marker))
                    elif projection:
                        axs[0].plot(x, y, color="r")
                        axs[1].plot(x, z, color="r")
                        axs[2].plot(z, y, color="r")
                    else:
                        if ((coords[2] - coords[3]) < slice_) and ((coords[2] + coords[3]) > slice_):
                            idx = ((x - coords[0])**2 + (y - coords[1])**2 < coords[3]**2 - (slice_ - coords[2])**2)
                            x = x[idx]
                            y = y[idx]
                            x = np.squeeze(x.reshape(1, -1))
                            y = np.squeeze(y.reshape(1, -1))
                            z = np.ones_like(x)*slice_
                            if z[0] == slice_:
                                axs[j].plot(x, y, color="darkgray")
                        else:
                            z = np.squeeze(z.reshape(1, -1))
    if not plot_3d:
        distance_from_borders = 0.0007
    # if plot_3d:
    #     layout = go.Layout(
    #         title='Wireframe Plot',
    #         scene=dict(
    #             xaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             ),
    #             yaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             ),
    #             zaxis=dict(
    #                 gridcolor='rgb(255, 255, 255)',
    #                 zerolinecolor='rgb(255, 255, 255)',
    #                 showbackground=True,
    #                 backgroundcolor='rgb(230, 230,230)'
    #             )
    #         ),
    #         showlegend=False,
    #     )
    #     fig = go.Figure(data=lines_plot, layout=layout)
    #     fig.show()
    if plot_traj:
        with open(traj_file) as f:
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
    if plot_3d:
        ax.set_xlim([0, max_lim])
        ax.set_ylim([0, max_lim])
        ax.set_zlim([0, max_lim])
    elif projection:
        for k in range(3):
            axs[k].set_xlim([0, max_lim])
            axs[k].set_ylim([0, max_lim])
    else:
        for k in range(len(z_slice)):
            axs[k].set_xlim([0, max_lim])
            axs[k].set_ylim([0, max_lim])
    plt.show()