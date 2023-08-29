import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from pathlib import Path
import glob

cur_path = os.getcwd()
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"

def get_dwi_array(dwi_path):
    # create an array with dwi values
    signal = []
    with open(dwi_path) as f:
        # read each line in the file
        # convert the line to a float and append it to the signal array
        [signal.append(float(line)) for line in f.readlines()]
    return np.array(signal)


def get_psge_data():
    # create an empty DataFrame with specified column names
    data_dwi = pd.DataFrame(columns=["x", "y", "z", "G", "Delta", "delta", "TE"])
    
    # create empty lists for each column
    x, y, z, G, Delta, delta, TE = [], [], [], [], [], [], []
    
    with open(scheme_file) as f:
        for line in f.readlines():
            if len(line.split(' ')) > 2:
                for e, element in enumerate(line.split(' ')):
                    if e == 0:
                        x.append(float(element))  # add x-coordinate value to the x list
                    elif e == 1:
                        y.append(float(element))  # add y-coordinate value to the y list
                    elif e == 2:
                        z.append(float(element))  # add z-coordinate value to the z list
                    elif e == 3:
                        G.append(float(element) * 1e-3)  # add G value (multiplied by 1e-3) to the G list
                    elif e == 4:
                        Delta.append(float(element))  # add Delta value to the Delta list
                    elif e == 5:
                        delta.append(float(element))  # add delta value to the delta list
                    elif e == 6:
                        TE.append(float(element[:-1]))  # add TE value (excluding the last character) to the TE list
    
    # assign the lists to respective columns in the DataFrame
    data_dwi["x"] = x
    data_dwi["y"] = y
    data_dwi["z"] = z
    data_dwi["G"] = G
    data_dwi["Delta"] = Delta
    data_dwi["delta"] = delta
    data_dwi["TE"] = TE
    
    # calculate a new column "b [ms/um²]" based on existing columns
    data_dwi["b [ms/um²]"] = pow(data_dwi["G"] * giro * data_dwi["delta"], 2) * \
                             (data_dwi["Delta"] - data_dwi["delta"] / 3) / 1000
    
    return data_dwi


def create_data(dwi_path, b0):
    # Get the DWI signal array from the given file path
    dwi_signal = get_dwi_array(dwi_path)
    
    # Get the psge data DataFrame
    data_psge = get_psge_data()
    
    # Add the DWI signal array as a new column named "DWI" in the data_psge DataFrame
    data_psge["DWI"] = list(dwi_signal)

    # Filter data based on x, y, and z coordinates greater than 0.0
    data_x = data_psge.loc[data_psge['x'] > 0.0]
    data_y = data_psge.loc[data_psge['y'] > 0.0]
    data_z = data_psge.loc[data_psge['z'] > 0.0]

    # Create a list of filtered dataframes
    datas = [data_x, data_y, data_z]

    for i in range(len(datas)):
        # Round b values to 1 decimal place
        datas[i]["b [ms/um²]"] = [round(n, 1) for n in list(datas[i]["b [ms/um²]"])]

        # Get the DWI(b0) value for the given b0
        Sb0 = list(datas[i].loc[datas[i]["b [ms/um²]"] == b0]["DWI"])[0]

        # Calculate log(Sb/So) for each DWI value in the dataframe
        signal = list(map(lambda Sb: np.log(Sb / Sb0), list(datas[i]["DWI"])))

        # Calculate adc [um²/ms] for each DWI value in the dataframe
        adc = list(map(lambda b, Sb: -np.log(Sb / Sb0) / (b - b0) if b != b0 else np.nan,
                       list(datas[i]["b [ms/um²]"]), list(datas[i]["DWI"])))

        # Add the calculated columns to the dataframe
        datas[i]["log(Sb/So)"] = signal
        datas[i]["adc [um²/ms]"] = adc

    # Concatenate the filtered dataframes into a single dataframe
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

def get_files_from_folder(folder_path):

    txt_files = glob.glob(os.path.join(folder_path, f"*.txt"))
    txt_files.sort()

    return txt_files