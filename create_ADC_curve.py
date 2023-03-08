# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"
icvf = 0.7
def get_dwi_array(dwi_path):
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
    data_dwi["b [ms/um²]"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2)*(data_dwi["Delta"]-data_dwi["delta"]/3)/1000

    return data_dwi

def create_data(dwi_path):
    dwi_signal = get_dwi_array(dwi_path)
    data_psge = get_psge_data()
    data_psge["DWI"] = list(dwi_signal)

    #x1 = list(data_psge["x"])[0]
    data_x = data_psge.loc[data_psge['x'] > 0.0]
    data_y = data_psge.loc[data_psge['y'] > 0.0]
    data_z = data_psge.loc[data_psge['z'] > 0.0]
    #data_1 = data_psge.loc[data_psge['x'] != x1]
    datas = [data_x,data_y,data_z]
    for i in range(len(datas)):
        b0 = list(datas[i]["b [ms/um²]"])[0]
        Sb0 = list(datas[i].loc[datas[i]["b [ms/um²]"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : np.log(Sb/Sb0), list(datas[i]["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(datas[i]["b [ms/um²]"]),list(datas[i]["DWI"])))
        datas[i]["log(Sb/So)"] = signal
        datas[i]["adc [um²/ms]"] = adc

    data_dwi = pd.concat(datas)

    return data_dwi

def get_adc(dwi_path):
    data = create_data(dwi_path)
    axes = ["x","y","z"]
    orientations = ["radial","radial","axial"]
    adcs = []
    for a in axes:
        ax_data = data.loc[data[a]>0]
        b1 = list(ax_data["b [ms/um²]"])[-1]
        ax_data = ax_data.loc[ax_data["b [ms/um²]"] == b1]
        adcs.append(list(ax_data["adc [um²/ms]"])[0])
    new_data = pd.DataFrame(columns = ["axis", "adc [um²/ms]"])
    new_data["orientations"] = orientations
    new_data["axis"] = axes
    new_data["adc [um²/ms]"] = adcs
    return new_data



def assemble_data(intra_active_path, intra_rest_path, extra_active_path, extra_rest_path):

    intra_active = get_adc(intra_active_path)
    intra_active["state"] = ["active"]*len(intra_active)
    intra_active["loc"] = ["intra"]*len(intra_active)
    extra_active = get_adc(extra_active_path)
    extra_active["state"] = ["active"]*len(extra_active)
    extra_active["loc"] = ["extra"]*len(extra_active)
    intra_rest = get_adc(intra_rest_path)
    intra_rest["state"] = ["rest"]*len(intra_rest)
    intra_rest["loc"] = ["intra"]*len(intra_rest)
    extra_rest = get_adc(extra_rest_path)
    extra_rest["state"] = ["rest"]*len(extra_rest)
    extra_rest["loc"] = ["extra"]*len(extra_rest)

    result = pd.concat([intra_active, intra_rest, extra_active, extra_rest])

    return result

def plot(result):
    fig1 = plt.figure(0)
    sns.catplot(x="orientations", y="adc [um²/ms]",
             hue="state", col = "loc", kind = "box",
             data=result)

    fig1 = plt.figure(1)
    sns.catplot(x="orientations", y="adc [um²/ms]",
             hue="state", col = "loc", 
             data=result)
    plt.show() 

def get_axons_array(file):
    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 3):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                if (len(spheres)!= 0):
                    axons.append(spheres)
                spheres = []
            

    return axons

def plot_radii(file):
    axons = get_axons_array(file)
    radiis =[] 
    area = 0
    for i in range(len(axons)):
        radius = [] 
        spheres = axons[i]
        for ii in range(len(spheres)):
            radius.append(spheres[ii][3])
        area += np.pi*spheres[ii][3]*spheres[ii][3]
        plt.plot(radius)
    print("approximate ICVF :"+str(area/(97.55*97.55)))
    plt.show()


def plot_tortuosity(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
        df = pd.DataFrame()
        df["x[um]"]=xs
        df["y[um]"]=ys
        df["index"] = [i]*len(xs)   
        dfs.append(df)
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="x[um]", y="y[um]")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()

def plot_tortuosity_angle(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        zs =[] 
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
            zs.append(spheres[ii][2]- spheres[ii-1][2])

        df = pd.DataFrame()
        df["index"] = [i]*len(xs)
        df["theta"]= list(map(lambda x,y:np.arctan(y/x)/np.pi if x!= 0 else np.pi/2, xs, ys))
        df["phi"]= list(map(lambda x,y,z:np.arctan(np.sqrt(x**2+y**2)/z)/np.pi if z!= 0 else np.pi/2, xs, ys, zs))  
        dfs.append(df)
        del df
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="theta", y="phi")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()
    

def plot_DWI(intra_rest_path,extra_rest_path, extra_active_path, intra_active_path):
    data_intra_rest = create_data(intra_rest_path)
    data_extra_active = create_data(extra_active_path)
    data_intra_active = create_data(intra_active_path)
    data_extra_rest = create_data(extra_rest_path)
    data_rest = data_intra_rest.add(data_extra_rest)
    data_active = data_intra_active.add(data_extra_active)
    ax = sns.lineplot(data=data_rest, x = "b [ms/um²]", y= "DWI")

    ax = sns.lineplot(data=data_active, x = "b [ms/um²]", y= "DWI")
    plt.show()





axons_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/axons_dist_50_0.3_gamma_distributed_axon_list.txt"

intra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/axons_dist_50_0.3_intra_active_DWI.txt"
extra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/axons_dist_50_0.3_extra_active_DWI.txt"
intra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/axons_dist_50_0.3_intra_rest_DWI.txt"
extra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/axons_dist_50_0.3_extra_rest_DWI.txt"               
data2 = assemble_data(intra_active_path, intra_rest_path, extra_active_path, extra_rest_path)

print(data2.groupby(["state"] ).mean())
sns.pointplot(data = data2.groupby(["state"]).mean().reset_index(), x = "adc [um²/ms]", y = "state")
plt.show()

plot_DWI(intra_rest_path,extra_rest_path, extra_active_path, intra_active_path)
