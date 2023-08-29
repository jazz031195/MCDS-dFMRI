import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
from read_outputs_functions import *

path = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/tortuous/icvf_0.50/"
folders = [f"{path}C_10_6/",f"{path}C_10_7/", f"{path}C_10_8/", f"{path}C_10_9/" ]

concentrations = [1e6, 1e7, 1e8, 1e9]
dfs = []
for e, folder in enumerate(folders):

    files = get_files_from_folder(folder)
    files = [file for file in files if "info" not in file]


    for file in files :

        if "T_1000" in file :
            df = get_adc(file, 1, 0.2)
            df["num_time_steps"] = [1000]*len(df)
            if "extra" in file:
                df["location"] = ["extra"]*len(df)
            elif "intra" in file:
                df["location"] = ["intra"]*len(df)
            df["concentration"] = [concentrations[e]]*len(df)
            dfs.append(df)
        if "T_5000" in file :
            df = get_adc(file, 1, 0.2)
            df["num_time_steps"] = [5000]*len(df)
            if "extra" in file:
                df["location"] = ["extra"]*len(df)
            elif "intra" in file:
                df["location"] = ["intra"]*len(df)
            df["concentration"] = [concentrations[e]]*len(df)
            dfs.append(df)
        elif "T_10000" in file :
            df = get_adc(file, 1, 0.2)
            df["num_time_steps"] = [10000]*len(df)
            if "extra" in file:
                df["location"] = ["extra"]*len(df)
            elif "intra" in file:
                df["location"] = ["intra"]*len(df)
            df["concentration"] = [concentrations[e]]*len(df)
            dfs.append(df)
        elif "T_30000" in file :
            df = get_adc(file, 1, 0.2)
            df["num_time_steps"] = [30000]*len(df)
            if "extra" in file:
                df["location"] = ["extra"]*len(df)
            elif "intra" in file:
                df["location"] = ["intra"]*len(df)
            df["concentration"] = [concentrations[e]]*len(df)
            dfs.append(df)


data = pd.concat(dfs)
print(data)
"""
g = sns.FacetGrid(data, col="orientations", row="location", margin_titles=True)
g.map_dataframe(sns.scatterplot, x="num_time_steps", y="adc [um²/ms]", hue="axis")
g.set_axis_labels("Num Time Steps", "ADC [um²/ms]")
g.add_legend()
plt.show()
"""

plt.figure(figsize=(8, 6))
sns.scatterplot(x="concentration", y="adc [um²/ms]", data=data.loc[data["location"]=="intra"], hue="num_time_steps", style="axis")
plt.xscale("log")  # Set x-axis to log scale
plt.title("ADC vs Concentration for Different Orientations : Intra")
plt.xlabel("Concentration (particles/mm3)")
plt.ylabel("ADC [um²/ms]")
plt.show()

plt.figure(figsize=(8, 6))
sns.scatterplot(x="concentration", y="adc [um²/ms]", data=data.loc[data["location"]=="extra"], hue="num_time_steps", style="axis")
plt.xscale("log")  # Set x-axis to log scale
plt.title("ADC vs Concentration for Different Orientations : Extra")
plt.xlabel("Concentration (particles/mm³)")
plt.ylabel("ADC [um²/ms]")
plt.show()