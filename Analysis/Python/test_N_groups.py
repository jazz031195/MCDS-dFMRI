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
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_21_dir.scheme"
icvf = 0.38
def get_dwi_array(dwi_path):
    # create array with dwi values
    signal = []
    with open(dwi_path) as f:
        [signal.append(float(line)) for line in f.readlines()]
    return np.array(signal)

def get_psge_data():
    data_dwi = pd.DataFrame(columns = ["x", "y","z","G","Delta","delta","TE"])
    x, y, z, G, Delta, delta, TE = [], [], [], [], [], [], []
    with open(scheme_file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 3:
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
    data_dwi["b [ms/um²]"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2) * (data_dwi["Delta"]-data_dwi["delta"]/3)/1000

    return data_dwi

def create_data(dwi_path):
    dwi_signal = get_dwi_array(dwi_path)
    data_psge = get_psge_data()
    data_psge["DWI"] = list(dwi_signal)
    nb_G = len(data_psge["G"].unique())
    nb_dir = int(len(dwi_signal) / nb_G)
    data_dwi = pd.DataFrame()
    for i in range(nb_dir):
        data_dir = data_psge.iloc[i*nb_G:(i+1)*nb_G]
        b0 = list(data_dir["b [ms/um²]"])[0]
        Sb0 = list(data_dir.loc[data_dir["b [ms/um²]"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : Sb/Sb0, list(data_dir["DWI"])))
        signal_log = list(map(lambda Sb : np.log(Sb/Sb0), list(data_dir["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(data_dir["b [ms/um²]"]),list(data_dir["DWI"])))
        data_dir["Sb/So"] = signal
        data_dir["log(Sb/So)"] = signal_log
        data_dir["adc [um²/ms]"] = adc
        data_dwi = pd.concat([data_dwi, data_dir])
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
                if t in T and n in N:
                    dwi_intra = create_data(DWI_folder / filename)
                    nb_G   = len(dwi_intra["G"].unique())
                    nb_dir = int(len(dwi_intra["x"].values) / nb_G)
                    for i in range(nb_G):
                        sb_so = []
                        for j in range(nb_dir):
                            sb_so.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So"])
                            b_lab = dwi_intra.iloc[nb_G*j + i, :]["b [ms/um²]"]
                        mean = np.mean(sb_so)
                        b_labels = dwi_intra["b [ms/um²]"].unique()
                        d = {'loc': "intra", 'N': n, 'T': t, 'Sb/So': mean, 'b [ms/um²]': b_lab}
                        df_avg_dwi = pd.DataFrame(d, index=[i])
                        df_dwi = pd.concat([df_dwi, df_avg_dwi])

    return df_dwi, df_crossings






branch = 'no_branching'
# combine_intra_extra_adc("neurons")
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/21_dir_benchmark/" + branch)

plot = False
T = ['15000'] #['1000', '5000', '10000', '15000']
N = ['10000', '15000', '30000', '55000']
df_dwi, _ = create_df(DWI_folder, T, N)
print(df_dwi)
df_dwi = df_dwi[df_dwi['b [ms/um²]'] > 0]
b_labels = df_dwi['b [ms/um²]'].unique()

fig, ax = plt.subplots(1, 1)
sns.violinplot(data=df_dwi, y='Sb/So', x='b [ms/um²]', hue='N', hue_order=N, ax=ax)
sns.move_legend(ax, "upper right")
ax.set_xticklabels([f'{float(blab):.3f}' for blab in b_labels])

couples = []
for b in df_dwi['b [ms/um²]'].unique():
    for N in df_dwi['N'].unique():
        couples.append((b, N))            

couples_end = []
for b in df_dwi['b [ms/um²]'].unique():
    for i in range(len(couples)):  
        for j in range(1, len(N)):
            if (i + j < len(couples)) and (couples[i][0] == b) and (couples[i+j][0] == b): 
                couples_end.append((couples[i], couples[i+j]))

statannot.add_stat_annotation(
    ax,
    data=df_dwi,
    y='Sb/So', x='b [ms/um²]',
    hue='N',
    hue_order=['10000', '15000', '30000', '55000'],
    box_pairs=couples_end,
    test="Mann-Whitney",
    text_format="star",
    loc="inside"
    )
ax.set_title(f"T = {T[0]}")
plt.show()

