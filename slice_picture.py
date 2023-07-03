import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm
import plotly.graph_objects as go
import plotly.colors as colors
from sklearn.neighbors import KNeighborsClassifier

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

def traj_data(name, zmin, zmax):
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

                indexes = [i for i,z in enumerate(zp) if z > zmin and z < zmax]

                xp = [xp[i] for i in indexes]
                yp = [yp[i] for i in indexes]
                zp = [zp[i] for i in indexes]

        return xp, yp, zp

def draw_cercles_with_trajectory(cylinder_array, size, file_name, zmin, zmax, swell = False):
    r = []
    fig, ax = plt.subplots(ncols = 2) 

    print(cylinder_array)

    xp, yp, zp = traj_data(file_name, zmin, zmax)

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

def main(name):

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
    
    draw_cercles(np.array(axons_slice), size, swell = False)


def main2(name, traj_name):
    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/{name}.txt"
    size = get_size (file)
    axons = get_spheres_array(file)
    axons_slice = []
    z = 0.01
    depth = 0.001
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
    
    draw_cercles_with_trajectory(np.array(axons_slice), size, traj_name, z-depth/2, z+depth/2, swell = False)

def split_list(lst, num_sublists):
    sublist_size = len(lst) // num_sublists
    remainder = len(lst) % num_sublists
    
    sublists = []
    start = 0
    
    for i in range(num_sublists):
        sublist_length = sublist_size + 1 if i < remainder else sublist_size
        sublists.append(lst[start:start+sublist_length])
        start += sublist_length
    
    return sublists

def main3(traj_name, zmin, zmax):


    colours = colors.qualitative.Plotly[:10]
    xs,ys,zs= traj_data(traj_name, zmin, zmax)

    # Create a scatter plot for the first DataFrame
    scatter = go.Scatter3d(
        x=xs,
        y=ys,
        z=zs,
        mode="markers",
        name=f"DataFrame ",
        marker=dict(
            size=1,  # Set the size of scatter points for DataFrame 1
            color=colours[0] # Set the color of scatter points for DataFrame 1
        )
    )

    # Create the layout
    layout = go.Layout(
        scene=dict(
            aspectmode="data"
        )
    )

    # Create the figure
    fig = go.Figure(data=scatter, layout=layout)


    # Show the figure
    fig.show()

def main4(name):
    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/{name}.txt"
    axons = get_spheres_array(file)
    scatters = []
    colours = colors.qualitative.Plotly[:10]
    for e, axon in enumerate(axons) :

        c = colours[e%10]
        # Create a scatter plot for the first DataFrame
        scatter = go.Scatter3d(
            x=[s[0]*1000 for s in axon],
            y=[s[1]*1000 for s in axon],
            z=[s[2]*1000 for s in axon],
            mode="markers",
            name=f"Axon {e}",
            marker=dict(
                sizemode="diameter",
                size=[s[3]*2000 for s in axon],  # Set the size of scatter points for DataFrame 1
                color=c, # Set the color of scatter points for DataFrame 1
                line=dict(
                    color="rgba(0, 0, 0, 0)",  # Set color to transparent (alpha=0)
                    width=0  # Set width to 0 to remove the contour
                )
            )
            
        )
        scatters.append(scatter)

    layout = go.Layout(
    scene=dict(
        xaxis=dict(title='X [um]'),
        yaxis=dict(title='Y [um]'),
        zaxis=dict(title='Z [um]')
    )
)

    # Create the figure
    fig = go.Figure(data=scatters, layout=layout)

    # Show the figure
    fig.show()

def main5 (traj_name_intra, traj_name_extra, zmin, zmax):

    xs,ys,zs= traj_data(traj_name_intra, zmin, zmax)

    # Create a scatter plot for the first DataFrame
    scatter_intra = go.Scatter3d(
        x=xs,
        y=ys,
        z=zs,
        mode="markers",
        name=f"DataFrame ",
        marker=dict(
            size=1,  # Set the size of scatter points for DataFrame 1
            color="green" # Set the color of scatter points for DataFrame 1
        )
    )
    xs,ys,zs= traj_data(traj_name_extra, zmin, zmax)

    # Create a scatter plot for the first DataFrame
    scatter_extra = go.Scatter3d(
        x=xs,
        y=ys,
        z=zs,
        mode="markers",
        name=f"DataFrame ",
        marker=dict(
            size=1,  # Set the size of scatter points for DataFrame 1
            color="red" # Set the color of scatter points for DataFrame 1
        )
    )

    # Create the layout
    layout = go.Layout(
        scene=dict(
            aspectmode="data"
        )
    )

    # Create the figure
    fig = go.Figure(data=[scatter_intra, scatter_extra], layout=layout)

    # Show the figure
    fig.show()

#main4( "icvf_0.6_swell_1.0_gamma_distributed_axon_list")
# main3("straight_try", 0,0.03)
main4("icvf_0.7_swell_1.0_gamma_distributed_axon_list")