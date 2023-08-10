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
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_21_dir.scheme"
# scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"
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
    data_dwi["b [um²/ms]"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2) * (data_dwi["Delta"]-data_dwi["delta"]/3)/1000

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
        b0 = list(data_dir["b [um²/ms]"])[0]
        Sb0 = list(data_dir.loc[data_dir["b [um²/ms]"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : Sb/Sb0, list(data_dir["DWI"])))
        signal_log = list(map(lambda Sb : np.log(Sb/Sb0), list(data_dir["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(data_dir["b [um²/ms]"]),list(data_dir["DWI"])))
        data_dir["Sb/So"] = signal
        data_dir["log(Sb/So)"] = signal_log
        data_dir["adc [um²/ms]"] = adc
        data_dwi = pd.concat([data_dwi, data_dir])
    return data_dwi


def create_df(DWI_folder):

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

            dwi_intra = create_data(DWI_folder / filename)
            nb_G   = len(dwi_intra["G"].unique())
            nb_dir = int(len(dwi_intra["x"].values) / nb_G)
            for i in range(nb_G):
                sb_so = []
                for j in range(nb_dir):
                    sb_so.append(dwi_intra.iloc[nb_G*j + i, :]["Sb/So"])
                    b_lab = dwi_intra.iloc[nb_G*j + i, :]["b [um²/ms]"]
                mean = np.mean(sb_so)
                b_labels = dwi_intra["b [um²/ms]"].unique()
                d = {'loc': "intra", 'N': n, 'T': t, 'Sb/So': mean, 'b [um²/ms]': b_lab}
                df_avg_dwi = pd.DataFrame(d, index=[i])
                df_dwi = pd.concat([df_dwi, df_avg_dwi])

    return df_dwi, df_crossings







# combine_intra_extra_adc("neurons")
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/21_dir/soma_only")
# DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/3_dir/soma_only")

plot = True
log  = False
df_dwi, df_crossings = create_df(DWI_folder)

T_labels = df_dwi['T'].unique()
N_labels = df_dwi['N'].unique()
b_labels = df_dwi["b [um²/ms]"].unique()
print(b_labels)
means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()

N_indices = np.argsort(N_labels.astype(float))
T_indices = np.argsort(T_labels.astype(float))
T_labels  = T_labels[T_indices]
N_labels  = N_labels[N_indices]
if log:
    y_lim_min = -5
    y_lim_max = 0.1
else:
    y_lim_min = 0.
    y_lim_max = 1.1

if plot:
    fig, ax = plt.subplots(1, 1)

for t_i, t in enumerate(T_labels):
    for n_i, n in enumerate(N_labels):
        signal_tmp = []
        err_tmp    = []
        for group, data in means.groupby(['N', 'T', 'b [um²/ms]']):
            N1 = group[0]
            T1 = group[1]
            b1 = group[2]
            S1 = data[0]

            if t==T1 and n==N1:
                if log:
                    signal_tmp.append(np.log(S1))
                else:
                    signal_tmp.append(S1)

        
        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2[0]
            N2 = group2[0]
            T2 = group2[1]
            if t==T2 and n==N2:
                err_tmp.append(err)

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            b_labels_shifted = [b_lab for b_lab in b_labels]
            if(len(signal_tmp) > 0):
                lines = ax.errorbar(b_labels_shifted, signal_tmp, yerr=err_tmp, label=f"N {n}, soma only", fmt='.')
                ax.set_xlabel('b [um²/ms]')
                if log:
                    ax.set_ylabel('ln(S/S0)')
                else:
                    ax.set_ylabel('S/S0')
                ax.legend(loc=1)
                ax.set_ylim([y_lim_min, y_lim_max])

# combine_intra_extra_adc("neurons")
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/21_dir/dendrites_only")

plot = True
df_dwi, df_crossings = create_df(DWI_folder)

T_labels = df_dwi['T'].unique()
N_labels = df_dwi['N'].unique()
b_labels = df_dwi["b [um²/ms]"].unique()
means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()

N_indices = np.argsort(N_labels.astype(float))
T_indices = np.argsort(T_labels.astype(float))
T_labels  = T_labels[T_indices]
N_labels  = N_labels[N_indices]

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

            if log:
                signal_tmp.append(np.log(S1))
            else:
                signal_tmp.append(S1)
        
        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2[0]
            N2 = group2[0]
            T2 = group2[1]
            if t==T2 and n==N2:
                err_tmp.append(err)

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            b_labels_shifted = [b_lab + 2*0.05 for b_lab in b_labels]
            if(len(signal_tmp) > 0):
                lines = ax.errorbar(b_labels_shifted, signal_tmp, yerr=err_tmp, label=f"N {n}, dendrites only", fmt='.')
                ax.legend(loc=1)
                ax.set_ylim([y_lim_min, y_lim_max])

# combine_intra_extra_adc("neurons")
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/21_dir/soma_dendrites")

plot = True
df_dwi, df_crossings = create_df(DWI_folder)

T_labels = df_dwi['T'].unique()
N_labels = df_dwi['N'].unique()
b_labels = df_dwi["b [um²/ms]"].unique()
means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()


N_indices = np.argsort(N_labels.astype(float))
T_indices = np.argsort(T_labels.astype(float))
T_labels  = T_labels[T_indices]
N_labels  = N_labels[N_indices]


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

            if log:
                signal_tmp.append(np.log(S1))
            else:
                signal_tmp.append(S1)

        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2[0]
            N2 = group2[0]
            T2 = group2[1]
            if t==T2 and n==N2:
                err_tmp.append(err)

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            b_labels_shifted = [b_lab + 4*0.05 for b_lab in b_labels]
            if(len(signal_tmp) > 0):
                lines = ax.errorbar(b_labels_shifted, signal_tmp, yerr=err_tmp, label=f"N {n}, soma & dendrites", fmt='.')
                ax.legend(loc=1)
                ax.set_ylim([y_lim_min, y_lim_max])

DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/21_dir/soma_dendrites_ex")

plot = True
df_dwi, df_crossings = create_df(DWI_folder)

T_labels = df_dwi['T'].unique()
N_labels = df_dwi['N'].unique()
b_labels = df_dwi["b [um²/ms]"].unique()
means = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].mean()
stds  = df_dwi.groupby(['N', 'T', 'b [um²/ms]'])['Sb/So'].std()

N_indices = np.argsort(N_labels.astype(float))
T_indices = np.argsort(T_labels.astype(float))
T_labels  = T_labels[T_indices]
N_labels  = N_labels[N_indices]

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

            if log:
                signal_tmp.append(np.log(S1))
            else:
                signal_tmp.append(S1)

        for group2, data2 in stds.groupby(['N', 'T', 'b [um²/ms]']):
            err = data2[0]
            N2 = group2[0]
            T2 = group2[1]
            if t==T2 and n==N2:
                err_tmp.append(err)

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            b_labels_shifted = [b_lab + 5*0.05 for b_lab in b_labels]
            if(len(signal_tmp) > 0):
                lines = ax.errorbar(b_labels_shifted, signal_tmp, yerr=err_tmp, label=f"N {n}, soma & dendrites ex", fmt='2', color='g')
                ax.legend(loc=1)
                ax.set_ylim([y_lim_min, y_lim_max])

# Analytical solutions & Mesh
if plot:
    # Analytical solutions
    G         = np.array([0, 0.015, 0.034, 0.048, 0.059, 0.107]) # in T/m
    Delta     = np.array([0.05] * G.size)  # in s
    delta     = np.array([0.0165] * G.size)# in s
    TE        = np.array([0.067] * G.size) # in s
    r_soma    = 10e-6 #m
    r_neurite = 0.5e-6 #m
    D0        = 2.5e-9 #m²/s
    gamma     = 2.6751525e8 #rad/(s*T)
    bb        = gamma**2 * G**2 * delta**2 * (Delta - delta/3) # rad² * s / m²
    print((D0/(gamma*G))**(1/3))

    nb_neurites     = 20
    l_neurite       = 240e-6 # m
    volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
    volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
    volume_neuron   = volume_neurites + volume_soma
    neurite_fraction= volume_neurites / volume_neuron
    soma_fraction   = volume_soma / volume_neuron
    print("{:e}".format((volume_neuron*1e9)))
    print("{:e}".format(1e10 * (volume_neuron*1e9)))

    soma_signal   = []
    soma_signal_neuman   = []
    neurites_signal = []
    both_signal     = []
    for i in range(bb.size):
        mlnS, mlnSneuman, mlnSnp, bardelta, b = my_murdaycotts(Delta[i], delta[i], r_soma, D0, bb[i])
        if log:
            soma_signal.append(-mlnS)
            soma_signal_neuman.append(-mlnSneuman)
        else:
            soma_signal.append(math.exp(-mlnS))
            soma_signal_neuman.append(math.exp(-mlnSneuman))

        # If b is 0, signal is 1 (no diffusion => no attenuation)
        if bb[i] == 0:
            if log:
                neurites_signal.append(0)
                both_signal.append(0)
            else:
                neurites_signal.append(1)
                both_signal.append(1)
        else:
            # Signal intra sticks
            if log:
                Ain = np.log(np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0)))
            else:
                Ain = np.sqrt(np.pi/(4 * bb[i] * D0)) * math.erf(np.sqrt(bb[i] * D0))
            neurites_signal.append(Ain)
            if log:
                both_signal.append(neurite_fraction * Ain + soma_fraction * -mlnS)
            else:
                both_signal.append(neurite_fraction * Ain + soma_fraction * math.exp(-mlnS))


    i = 0
    for t_i, t in enumerate(T_labels):

        if np.sum(np.isnan(err_tmp)) == 0 and plot:
            ax2 = ax.twinx()
            # Replace the NaN corresponding to b=0 to 1
            ax2.errorbar([b_lab + 0.05 for b_lab in b_labels], soma_signal, 
                            yerr=[0], label=f"Soma (analytic)", fmt='*')
            ax2.errorbar([b_lab + 0.05 for b_lab in b_labels], soma_signal_neuman, 
                            yerr=[0], label=f"Soma (analytic, Neuman)", fmt='o', color='blue')
            ax2.errorbar([b_lab + 3*0.05 for b_lab in b_labels], neurites_signal, 
                            yerr=[0], label=f"Neurites (analytic)", fmt='*')
            ax2.errorbar([b_lab + 6*0.05 for b_lab in b_labels], both_signal, 
                            yerr=[0], label=f"Neurites & soma (analytic)", fmt='*')
            ax2.plot([0, 10], [0, -4], label="D = 2.5 [ms/um²]")
            ax2.legend(loc=3)
            ax2.set_yticklabels([])
            ax2.set_yticks([])
            ax2.set_ylim([y_lim_min, y_lim_max])
            step_length = np.sqrt(6 * D0 * TE[0] / int(t))
            ax2.set_title(f"T = {T_labels[i]}, step length = {step_length*1e6:.3f} um")
            i = i + 1
if log:
    fig.suptitle('ln(S/S0) average over 21 directions, average over 5 rep', y=0.95)
else:
    fig.suptitle('S/S0 average over 21 directions, average over 5 rep', y=0.95)
plt.show()
