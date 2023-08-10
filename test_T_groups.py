# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
import pandas as pd
import seaborn as sns
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings("ignore")
sys.path.insert(1, '/home/localadmin/Documents/analytical_formula/')
from my_murdaycotts import my_murdaycotts
import math
import statannot

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_3_dir_5_b.scheme"
icvf = 0.38
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
    data_dwi["b [um²/ms]"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2)*(data_dwi["Delta"]-data_dwi["delta"]/3)/1000

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
        b0 = list(datas[i]["b [um²/ms]"])[0]
        Sb0 = list(datas[i].loc[datas[i]["b [um²/ms]"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : Sb/Sb0, list(datas[i]["DWI"])))
        signal_log = list(map(lambda Sb : np.log(Sb/Sb0), list(datas[i]["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(datas[i]["b [um²/ms]"]),list(datas[i]["DWI"])))
        datas[i]["Sb/So"] = signal
        datas[i]["log(Sb/So)"] = signal_log
        datas[i]["adc [um²/ms]"] = adc

    data_dwi = pd.concat(datas)

    return data_dwi

def create_df(DWI_folder, T, N):

    df_dwi = pd.DataFrame()
    df_crossings = pd.DataFrame()

    # Iterate through the files in the folder
    for filename in os.listdir(DWI_folder):
        # Check if the filename contains "_rep_" and "DWI"
        if "DWI" in filename:
            # Split the string by dots ('.') to separate the filename and extension
            filename_without_ext, extension = filename.split('.')
            # Split the filename by underscores
            parts = filename_without_ext.split('_')
            parts = parts[:-1]
            # Add the desired suffix to the last element of the parts list
            parts[-1] = f"{parts[-1]}_simulation_info"
            # Join the parts back together using underscores
            filename_simu_info = '_'.join(parts) + '.' + extension

            # Initialize a variable to store the number of particles eliminated due to crossings
            num_particles_crossings = None

            # Open the file in read mode
            with open(DWI_folder / filename_simu_info, 'r') as file:
                # Read the file line by line
                for line in file:
                    if "Number of particles:" in line:
                        # :-1 to remove the /n at the end
                        n = line.split('-')[-1][:-1]
                    if "Number of steps:" in line:
                        t = line.split('-')[-1][:-1]
                    # Check if the line contains the relevant information
                    if 'Number of particles eliminated due crossings' in line:
                        # Split the line to get the number of particles as the last element
                        num_particles_crossings = int(line.split()[-1])
                        # Break the loop, as we have found the information we need
                        break
                d = {'nb_crossings': [num_particles_crossings], 'N': [n], 'T': [t]}
                df_avg_crossings = pd.DataFrame(d)
                df_crossings = pd.concat([df_crossings, df_avg_crossings])
                
                if t in T and n == N:
                    dwi_intra = create_data(DWI_folder / filename)
                    data_x = dwi_intra.loc[(dwi_intra['x'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
                    data_y = dwi_intra.loc[(dwi_intra['y'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
                    data_z = dwi_intra.loc[(dwi_intra['z'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
                    datas= [data_x, data_y, data_z]
                    mean = list(map(lambda x,y,z : np.mean([x,y,z]), list(datas[0]["Sb/So"]),list(datas[1]["Sb/So"]),list(datas[2]["Sb/So"])))
                    b_labels = data_x["b [um²/ms]"].unique()
                    d = {'loc': "intra", 'N': n, 'T': t, 'Sb/So': mean, 'b [um²/ms]': b_labels}
                    df_avg_dwi = pd.DataFrame(d)
                    df_dwi = pd.concat([df_dwi, df_avg_dwi])

    return df_dwi, df_crossings







# combine_intra_extra_adc("neurons")
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/3_dir/branching")

plot = False
T = ['1000', '5000', '10000', '15000']
N = str(10000)
df_dwi, _ = create_df(DWI_folder, T, N)
df_dwi = df_dwi[df_dwi['b [um²/ms]'] > 0]
b_labels = df_dwi['b [um²/ms]'].unique()

fig, ax = plt.subplots(1, 1)
sns.violinplot(data=df_dwi, y='Sb/So', x='b [um²/ms]', hue='T', hue_order=['1000', '5000', '10000', '15000'], ax=ax)
sns.move_legend(ax, "upper right")
ax.set_xticklabels([f'{float(blab):.3f}' for blab in b_labels])

couples = []
for b in df_dwi['b [um²/ms]'].unique():
    for T in df_dwi['T'].unique():
        couples.append((b, T))            

couples_end = []
for b in df_dwi['b [um²/ms]'].unique():
    for i in range(len(couples)):  
        if (i + 1 < len(couples)) and (couples[i][0] == b) and (couples[i+1][0] == b): 
            couples_end.append((couples[i], couples[i+1]))
        if (i + 2 < len(couples)) and (couples[i][0] == b) and (couples[i+2][0] == b): 
            couples_end.append((couples[i], couples[i+2]))
statannot.add_stat_annotation(
    ax,
    data=df_dwi,
    y='Sb/So', x='b [um²/ms]',
    hue='T',
    hue_order=['1000', '5000', '10000', '15000'],
    box_pairs=couples_end,
    test="Mann-Whitney",
    text_format="star",
    loc="inside"
    )

plt.show()

