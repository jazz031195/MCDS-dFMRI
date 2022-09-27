import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation
from scipy.linalg import norm

def data_for_cylinder(p0,p1,R):
    #vector in direction of axis
    v = p1 - p0
    #find magnitude of vector
    mag = norm(v)
    #unit vector in direction of axis
    v = v / mag
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0])
    if (v == not_v).all():
        not_v = np.array([0, 1, 0])
    #make vector perpendicular to v
    n1 = np.cross(v, not_v)
    #normalize n1
    n1 /= norm(n1)
    #make unit vector perpendicular to v and n1
    n2 = np.cross(v, n1)
    #surface ranges over t from 0 to length of axis and 0 to 2*pi
    t = np.linspace(0, mag, 100)
    theta = np.linspace(0, 2 * np.pi, 100)
    #use meshgrid to make 2d arrays
    t, theta = np.meshgrid(t, theta)
    #generate coordinates for surface
    X, Y, Z = [p0[i] + v[i] * t + R * np.sin(theta) * n1[i] + R * np.cos(theta) * n2[i] for i in [0, 1, 2]]
    return X, Y, Z

def get_T_N (cur_path):
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
    return T, N

def get_nbr_cylinders (cur_path):
    # find configuration file for N and T values
    conf_file_path = cur_path + "/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"
    nbr_cylinders = 0
    with open(conf_file_path) as f:
        for line in f.readlines():
            if line.split(' ')[0] == "num_cylinders":
                nbr_cylinders = int(line.split(' ')[1][:-1])
    return nbr_cylinders

def get_trajectory_list(cur_path):
    # find path with trajectories
    file_name = "cylinder_gamma_packing_test_0.traj.txt"
    traj_path = cur_path +"/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/"+ str(file_name)
    # create array with trajectory values
    contents = []
    with open(traj_path) as f:
        [contents.append(float(line)) for line in f.readlines()]
    return contents 

def get_trajectory_array(T,N, cur_path):
    contents = get_trajectory_list(cur_path)
    # initialise trajectory array
    traj_array = np.zeros((T+1, N, 3))
    # fill array with values of contents
    for i in range(len(contents)):
        traj_array[int(i/3)%(T+1)][int(int(i/3)/(N+1))][i%3] = contents[i]
    return traj_array

def get_cylinder_array (cur_path):
    nbr_cylinders = get_nbr_cylinders (cur_path)
    # find configuration file for N and T values
    conf_file_path = cur_path + "/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/cylinder_gamma_packing_test_gamma_distributed_cylinder_list.txt"
    cylinder_array = np.zeros((nbr_cylinders,7))
    with open(conf_file_path) as f:
        e = 0
        for line in f.readlines():
            if e< nbr_cylinders:
                if len(line.split(' ')) > 2:
                    cylinder_array[e] = np.array([float(i)*0.001 for i in line.split(' ')])
                    e = e+1
    return cylinder_array

def get_dwi_array(cur_path):
    # find path with trajectories
    file_name = "cylinder_gamma_packing_test_DWI.txt"
    traj_path = cur_path +"/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/"+ str(file_name)
    # create array with trajectory values
    signal = []
    with open(traj_path) as f:
        [signal.append(float(line)) for line in f.readlines()]
    return np.array(signal)


def main(Nmax = 5, Cmax= 5, plot_dwi = True):
    # find current path
    cur_path = os.getcwd()
    # find number of time steps (T) and number of walkers (N)
    T, N = get_T_N (cur_path)
    # find trajectory array
    traj_array = get_trajectory_array(T,N, cur_path)
    # find cylinder positions 
    cylinder_array = get_cylinder_array (cur_path)
    # animation 
    def update_graph(num):
        for n in range(Nmax):
            data = traj_array[num,n,:]
            graph[n]._offsets3d = (pd.Series(data[0]), pd.Series(data[1]), pd.Series(data[2]))
        title.set_text('3D Test, time={}'.format(num))

    fig1 = plt.figure(1)
    ax = fig1.add_subplot(111, projection='3d')
    title = ax.set_title('3D Test')

    data=traj_array[0]
    graph = [ ax.scatter(pd.Series(data[n][0]), pd.Series(data[n][1]), pd.Series(data[n][2])) for n in range(Nmax) ]

    ani = matplotlib.animation.FuncAnimation(fig1, update_graph, T, 
                                interval=5, blit=False)
    
    for c in range(Cmax):
        cylinder = cylinder_array[c]
        Xc,Yc,Zc = data_for_cylinder(np.array([cylinder[0],cylinder[1],cylinder[2]]), np.array([cylinder[3],cylinder[4],cylinder[5]*0.1]), R = cylinder[6])
        ax.plot_surface(Xc, Yc, Zc)
    

    if plot_dwi:
        dwi_signal = get_dwi_array(cur_path)
        fig2 = plt.figure(2)
        x = np.linspace(0,len(dwi_signal)-1, len(dwi_signal))
        plt.plot(x, dwi_signal)
        plt.title('DWI Signal')

    plt.show()
main(Nmax= 5, Cmax=5)

