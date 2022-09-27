import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation

def main(Nmax):

    # find current path
    cur_path = os.getcwd()
    # find path with trajectories
    file_name = "cylinder_gamma_packing_test_0.traj.txt"
    traj_path = cur_path +"/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/"+ str(file_name)
    # create array with trajectory values
    contents = []
    with open(traj_path) as f:
        [contents.append(float(line)) for line in f.readlines()]
    # find configuration file for N and T values
    conf_file_path = cur_path + "/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"
    N = 0
    T = 0
    with open(conf_file_path) as f:
        for line in f.readlines():
            if line.split(' ')[0] == "N":
                N = int(line.split(' ')[1][:-1])
            elif line.split(' ')[0] == "T":
                T = int(line.split(' ')[1][:-1])
    # initialise trajectory array
    traj_array = np.zeros((T+1, N, 3))
    # fill array with values of contents
    for i in range(len(contents)):
        traj_array[int(i/3)%(T+1)][int(int(i/3)/(N+1))][i%3] = contents[i]

    # animation 
    def update_graph(num):
        for n in range(Nmax):
            data = traj_array[num,n,:]
            graph[n]._offsets3d = (pd.Series(data[0]), pd.Series(data[1]), pd.Series(data[2]))
        title.set_text('3D Test, time={}'.format(num))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    title = ax.set_title('3D Test')

    data=traj_array[0]
    graph = [ ax.scatter(pd.Series(data[n][0]), pd.Series(data[n][1]), pd.Series(data[n][2])) for n in range(Nmax) ]

    ani = matplotlib.animation.FuncAnimation(fig, update_graph, T, 
                                interval=5, blit=False)

    plt.show()

main(25)