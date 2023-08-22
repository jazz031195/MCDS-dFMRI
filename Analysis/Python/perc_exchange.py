import re
from pathlib import Path
import os 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

branching = 'branching'
DWI_folder = Path("/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/exchange")
df = pd.DataFrame(columns = ["Soma begin", "Soma end", "Dendrites begin","Dendrites end", "neuron"], dtype="object")
for neuron in os.listdir(DWI_folder):
    for file in os.listdir(DWI_folder / neuron / "perc_exchange"):
        f = open(DWI_folder / neuron / f"perc_exchange/{file}", "r")
        numbers = re.findall(r'\d+', f.read())
        list_ = numbers[:4]
        list_.append(neuron)
        df.loc[len(df)] = list_


df["Soma begin"]      = df["Soma begin"].astype(int)
df["Soma end"]        = df["Soma end"].astype(int)
df["Dendrites begin"] = df["Dendrites begin"].astype(int)
df["Dendrites end"]   = df["Dendrites end"].astype(int)


sum = df.groupby(['neuron']).sum() / 5 

MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


perc  = sum / 50000 * 100
means = perc.mean().values
stds  = perc.std().values
plt.errorbar(["Soma begin", "Soma end", "Dendrites begin", "Dendrites end"], perc.mean().values, yerr=perc.std().values, fmt='.')
plt.ylabel('% walker')
plt.show()



ind = np.arange(len(means) / 2)  # the x locations for the groups
width = 0.35  # the width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, [means[0], means[2]], width, yerr=[stds[0], stds[2]],
                label='Begin')
rects2 = ax.bar(ind + width/2, [means[1], means[3]], width, yerr=[stds[1], stds[3]],
                label='End')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('% water molecules')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
ax.legend()
plt.show()


nb_neurites     = 20
r_soma    = 10e-6 #m
r_neurite = 0.5e-6 #m
if branching == 'branching':
    nb_branching = 3
    l_neurite       = 80e-6 #m
    volume_neurites = nb_neurites * (2**nb_branching + 1) * np.pi*r_neurite**2*l_neurite # in m³
else:
    l_neurite       = 240e-6 # m
    volume_neurites = nb_neurites * np.pi*r_neurite**2*l_neurite # in m³
volume_soma     = 4/3 * np.pi * r_soma**3 # in m³
volume_soma     = volume_soma * 1e18
volume_neurites = volume_neurites * 1e18
volume_neurites = 8782.71
volume_neuron   = volume_neurites + volume_soma

means = sum.mean().values
stds  = sum.std().values
print(means, volume_soma, volume_neurites)
means[[0, 1]] = means[[0, 1]] / volume_soma
means[[2, 3]] = means[[2, 3]] / volume_neurites
stds[[0, 1]]  = stds[[0, 1]]  / volume_soma
stds[[2, 3]]  = stds[[2, 3]]  / volume_neurites
# 4.18879e-06 8.78271e-06
print(means, volume_soma, volume_neurites)

fig, ax = plt.subplots()
rects1 = ax.bar(ind - width/2, [means[0], means[2]], width, yerr=[stds[0], stds[2]],
                label='Begin')
rects2 = ax.bar(ind + width/2, [means[1], means[3]], width, yerr=[stds[1], stds[3]],
                label='End')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Water molecules density [part/um³]')
ax.set_xticks(ind)
ax.set_xticklabels(('Soma', 'Dendrites'))
ax.legend()
plt.show()