# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

cur_path = os.getcwd()
N_5_10_4_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_5_10_4/dyn_cylinder_gamma_rest_extra_DWI.txt"
N_10_4_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_10_4/dyn_cylinder_gamma_rest_extra_DWI.txt"
N_5_10_5_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_5_10_5/dyn_cylinder_gamma_rest_extra_DWI.txt"
N_10_5_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_10_5/dyn_cylinder_gamma_rest_extra_DWI.txt"
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
scheme_file = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"
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
    data_dwi["b"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2)*(data_dwi["Delta"]-data_dwi["delta"]/3)

    return data_dwi

def create_data(dwi_path):
    dwi_signal = get_dwi_array(dwi_path)
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

    data_dwi = pd.concat([data_1, data_2])

    return data_dwi
    
def plot_DWI(data, path):
    fig1 = plt.figure(1)
    sns.lineplot(x="b", y="signal",
             hue="x", style="y",
             data=data)
    plt.title('DWI Signal (log(Sb/So))')
    path_ = path.split("/")
    path_ = path_[:-1]
    path_ = "/".join(path_)
    name = path_ + "/DWI.png"
    print(name)
    fig1.savefig(name)

    plt.show() 
def get_ADC_value(data_dwi):
    x0 = list(data_dwi["x"])[0]
    b_200 = list(data_dwi["b"])[1]
    b_1000 = list(data_dwi["b"])[-1]
    data_dwi0 = data_dwi.loc[data_dwi.x == x0]
    signalb200 = list(data_dwi0.loc[data_dwi0.b == b_200]["signal"])[0]
    signalb1000 = list(data_dwi0.loc[data_dwi0.b == b_1000]["signal"])[0]
    return (signalb200-signalb1000)/800

def main():
    adcs = []
    Ns = [1e4, 5e4, 1e5, 5e5]
    for path in [N_10_4_file,N_5_10_4_file, N_10_5_file, N_5_10_5_file]:
        data_dwi = create_data(path)
        plot_DWI(data_dwi, path)
        adcs.append(get_ADC_value(data_dwi))
    print(adcs)
    fig2 = plt.figure(2)
    plt.plot(Ns, adcs)
    path_ = path.split("/")
    path_ = path_[:-2]
    path_ = "/".join(path_)
    name = path_ + "/ADC_sanity_check.png"
    plt.title("ADC [mm²/s] wrt number of walkers (diffusion coefficient = 2.5e-3 mm²/s)")
    fig2.savefig(name)

    
main()





