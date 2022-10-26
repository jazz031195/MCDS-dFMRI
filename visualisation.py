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
file_name = "dyn_cylinder_gamma_rest_extra_DWI.txt"
dwi_path = cur_path +"/MCDC_Simulator_public-master/instructions/demos/output/"+ str(file_name)
gamma_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dyn_cylinder_gamma_packing_test_gamma_distributed_dyn_cylinder_list.txt"
file_name = "dyn_cylinder_gamma_packing_test_0.traj.txt"
gamma_traj_path = cur_path +"/MCDC_Simulator_public-master/instructions/demos/output/"+ str(file_name)
conf_file_path = cur_path + "/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedCylinders.conf"
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
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
    nbr_cylinders = 0
    with open(conf_file_path) as f:
        for line in f.readlines():
            if line.split(' ')[0] == "num_cylinders":
                nbr_cylinders = int(line.split(' ')[1][:-1])
    return nbr_cylinders

def get_trajectory_list(cur_path):
    # create array with trajectory values
    contents = []
    with open(gamma_traj_path) as f:
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
    cylinder_array = np.zeros((nbr_cylinders,7))
    with open(gamma_file_path) as f:
        e = 0
        for line in f.readlines():
            if e< nbr_cylinders:
                if len(line.split(' ')) > 2:
                    cylinder_array[e] = np.array([float(i)*0.001 for i in line.split(' ')])
                    e = e+1
    return cylinder_array

def get_dwi_array(cur_path):
    # create array with dwi values
    signal = []
    with open(dwi_path) as f:
        [signal.append(float(line)) for line in f.readlines()]
    return np.array(signal)

def get_psge_data():
    data_dwi = pd.DataFrame(columns = ["x", "y","z","G","Delta","delta","TE"])
    x, y, z,G,Delta,delta,TE = [],[],[],[],[],[],[]
    with open(scheme_file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 2:
                print(line.split(' '))
                for e, element in enumerate(line.split(' ')):
                    if e == 0:
                        x.append(float(element))
                    elif e == 1:
                        y.append(float(element))
                    elif e == 2:
                        z.append(float(element))
                    elif e == 3:
                        G.append(float(element)*1e-3)
                    elif e == 4:
                        Delta.append(float(element))
                    elif e == 5:
                        delta.append(float(element))
                    elif e == 6:
                        TE.append(float(element[:-1]))
    data_dwi["x"] = x
    data_dwi["y"] = y
    data_dwi["z"] = z
    data_dwi["G"] = G
    data_dwi["Delta"] = Delta
    data_dwi["delta"] = delta
    data_dwi["TE"] = TE
    data_dwi["b"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2)*(data_dwi["Delta"]-data_dwi["delta"]/3)

    return data_dwi





def animation(Nmax = 5, Cmax= 5):
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
        

def plot_adc():

    dwi_signal = get_dwi_array(cur_path)
    data_psge = get_psge_data()
    data_psge["DWI"] = list(dwi_signal)

    x1 = list(data_psge["x"])[0]
    data_1 = data_psge.loc[data_psge['x'] == x1]
    data_2 = data_psge.loc[data_psge['x'] != x1]

    Sb0_1 = list(data_1.loc[data_1["b"]== 0]["DWI"])[0]
    signal = list(map(lambda Sb : np.log(Sb/Sb0_1), list(data_1["DWI"])))
    data_1["signal"] = signal
    Sb0_2 = list(data_2.loc[data_2["b"]== 0]["DWI"])[0]
    signal = list(map(lambda Sb : np.log(Sb/Sb0_2), list(data_2["DWI"])))
    data_2["signal"] = signal

    data = pd.concat([data_1, data_2])
    fig2 = plt.figure(2)
    sns.lineplot(x="b", y="signal",
             hue="x", style="y",
             data=data)
    plt.title('DWI Signal [mm²/s]')
    fig2.savefig("ADC.png")

    plt.show()

plot_adc()


