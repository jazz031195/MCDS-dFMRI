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

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"
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

def get_adc(dwi_path, b, b0):
    data = create_data(dwi_path, b0)
    axes = ["x","y","z"]
    orientations = ["radial","radial","axial"]
    adcs = []
    for a in axes:
        ax_data = data.loc[data[a]>0]
        ax_data = ax_data.loc[ax_data["b [ms/um²]"] == b]
        adcs.append(list(ax_data["adc [um²/ms]"])[0])
    new_data = pd.DataFrame(columns = ["axis", "adc [um²/ms]"])
    new_data["orientations"] = orientations
    new_data["axis"] = axes
    new_data["adc [um²/ms]"] = adcs
    return new_data

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
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/")

plot = True
T = ['1000', '5000', '10000', '15000']
N = str(5000)
df_dwi, _ = create_df(DWI_folder, T, N)
T_labels = T
N_labels = [str(N)]
b_labels = df_dwi["b [um²/ms]"].unique()
means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()


if plot:
    fig, ax = plt.subplots(1, 1)

for t_i, t in enumerate(T_labels):
    for n_i, n in enumerate(N_labels):
        signal_tmp = []
        err_tmp    = []
        nb_crossings_tmp = []
        for group, data in means.groupby(['N', 'T', 'b [um²/ms]']):
            N1 = group[0]
            T1 = group[1]
            b1 = group[2]
            S1 = data[0]

            if t==T1 and n==N1:
                signal_tmp.append(S1)
        
        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2[0]
            N2 = group2[0]
            T2 = group2[1]
            if t==T2 and n==N2:
                err_tmp.append(err)


        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            lines = ax.errorbar(b_labels, signal_tmp, yerr=err_tmp, label=f"N {n}, T {t}", fmt='--')
            ax.set_xlabel('b values')
            ax.set_ylabel('S/S0')
            ax.set_ylim([0, 1.1])
            ax.legend(loc=1)

# Analytical solutions & Mesh
if plot: 
    G         = np.array([0, 0.015, 0.034, 0.048, 0.059]) # in T/m
    Delta     = np.array([0.05] * G.size)  # in s
    delta     = np.array([0.0165] * G.size)# in s
    TE        = np.array([0.067] * G.size) # in s
    r_soma    = 10e-6 #m
    r_neurite = 0.5e-6 #m
    D0        = 2.5e-9 #m²/s
    gamma     = 2.6751525e8 #rad/(s*T)
    bb        = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s / m²

    nb_neurites     = 20
    l_neurite       = 245 # m
    volume_neurites = nb_neurites * np.pi*(r_neurite*1e6)**2*l_neurite # in m³
    volume_soma     = 4/3 * np.pi * (r_soma*1e6)**3 # in m³
    volume_neuron   = volume_neurites + volume_soma
    neurite_fraction= volume_neurites / volume_neuron
    soma_fraction   = volume_soma / volume_neuron


    soma_signal     = []
    neurites_signal = []
    both_signal     = []
    for i in range(bb.size):
        mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(Delta[i], delta[i], r_soma, D0, bb[i])
        soma_signal.append(math.exp(-mlnS))

        # Signal intra sticks
        Ain = np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0))
        neurites_signal.append(Ain)
        both_signal.append(neurite_fraction * Ain + soma_fraction * math.exp(-mlnS))


    mesh_folders = ['soma', 'nb_20_1_subbranch_245_span']
    ax2 = ax.twinx()
    ax2.set_ylim([0, 1.1])
    # Replace the NaN corresponding to b=0 to 1
    neurites_signal[0] = both_signal[0] = 1

    ax2.errorbar(b_labels, soma_signal, 
                    yerr=[0], label=f"Soma (analytic)", fmt='-*')
    ax2.errorbar(b_labels, neurites_signal, 
                    yerr=[0], label=f"Neurites (analytic)", fmt='-*')
    ax2.errorbar(b_labels, both_signal, 
                    yerr=[0], label=f"Neurites & soma (analytic)", fmt='-*')
    for mesh in mesh_folders:
        DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/mesh") / mesh

        df_dwi, _ = create_df(DWI_folder, T, N)
        
        means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
        stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()
        
        signal_tmp = []
        err_tmp    = []
        for group, data in means.groupby(['N', 'T', 'b [um²/ms]']):
            N1 = group[0]
            T1 = group[1]
            b1 = group[2]
            S1 = data.values[0]
            if T1 == '5000':
                signal_tmp.append(S1)
        
        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2.values[0]
            T2 = group2[1]
            if T2 == '5000':
                err_tmp.append(err)

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            if signal_tmp:
                ax2.errorbar(b_labels, signal_tmp, yerr=err_tmp, label=f"N {n} {mesh} (mesh)", fmt='-o')
                
                

ax2.legend(loc=3)
ax2.set_yticklabels([])
ax2.set_yticks([])

fig.suptitle('S/S0 average over x, y, z direction, average over 10 rep', y=0.95)
plt.show()

