# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

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
    data_dwi["b"] = pow(data_dwi["G"]*giro*data_dwi["delta"],2)*(data_dwi["Delta"]-data_dwi["delta"]/3)

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
        Sb0 = list(datas[i].loc[datas[i]["b"]== 0]["DWI"])[0]
        signal = list(map(lambda Sb : np.log(Sb/Sb0), list(datas[i]["DWI"])))
        datas[i]["log(Sb/So)"] = signal

    data_dwi = pd.concat(datas)

    return data_dwi
    
def plot_DWI(data, path):
    fig1 = plt.figure(1)
    sns.lineplot(x="b", y="log(Sb/So)",
             hue="x", style="y",size = "z",
             data=data)
    plt.title('DWI Signal (log(Sb/So))')
    path_ = path.split("/")
    path_ = path_[:-1]
    path_ = "/".join(path_)
    name = path_ + "/DWI.png"

    fig1.savefig(name)

    plt.show() 

def get_ADC_value(data_dwi, orientation = "all"):
    xs= [1.0,0,0]   
    ys = [0,1.0,0]
    zs = [0,0,1.0]
    adcs =[] 
    for i in range(3):  
        x0 = xs[i] 
        y0 = ys[i] 
        z0 = zs[i] 
        b_200 = list(data_dwi["b"])[1]
        b_1000 = list(data_dwi["b"])[-1]
        data_dwi0 = data_dwi.loc[data_dwi.x == x0]
        data_dwi0 = data_dwi0.loc[data_dwi0.y == y0]
        data_dwi0 = data_dwi0.loc[data_dwi0.z == z0]
        signalb200 = list(data_dwi0.loc[data_dwi0.b == b_200]["log(Sb/So)"])[0]
        signalb1000 = list(data_dwi0.loc[data_dwi0.b == b_1000]["log(Sb/So)"])[0]
        adcs.append((signalb200-signalb1000)/(b_1000-b_200))
    if orientation == "all":
        return np.sum(adcs)/3
    if orientation == "separate":
        return adcs[0], adcs[1], adcs[2]


def get_total_dwi_array(path):
    
    data_dwi = create_data(path)
    data_x = data_dwi.loc[data_dwi['x'] > 0.0].sort_values(by = ["b"],ascending=True)
    data_y = data_dwi.loc[data_dwi['y'] > 0.0].sort_values(by = ["b"],ascending=True)
    data_z = data_dwi.loc[data_dwi['z'] > 0.0].sort_values(by = ["b"],ascending=True)
    bs = list(data_z["b"])
    datas = [data_x,data_y,data_z]

    #x0 = list(dict.fromkeys(list(data_dwi["x"])))[0]
    #data_dwi0 = list(data_dwi.loc[data_dwi.x == x0].sort_values(by = ["b"],ascending=True)["DWI"])
    #x1 = list(dict.fromkeys(list(data_dwi["x"])))[1]
    #data_dwi1 = list(["DWI"])
    #bs = list(data_dwi.loc[data_dwi.x == x1].sort_values(by = ["b"])["b"])
    dwi_tot = list(map(lambda x,y,z : pow(x*y*z,1/3), list(datas[0]["DWI"]),list(datas[1]["DWI"]),list(datas[2]["DWI"])))
    dwi_tot = list(map(lambda x : np.log(x/dwi_tot[0]), dwi_tot))
    return bs, dwi_tot



def plot_adc_sanity_check(state, loc):

    orientations = []
    adcs = []
    confs_ =[]  
    data_adc = pd.DataFrame()
    Ns_ =[] 
    Ns = [1e4, 5e4, 1e5, 5e5, 1e6]
    
    for conf in ["conf1","conf2","conf3"] :
        N_5_10_4_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/N_5_10_4/"+str(conf)+"/dyn_cylinder_DWI.txt"
        N_10_4_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/N_10_4/"+str(conf)+"/dyn_cylinder_DWI.txt"
        N_5_10_5_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/N_5_10_5/"+str(conf)+"/dyn_cylinder_DWI.txt"
        N_10_5_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/N_10_5/"+str(conf)+"/dyn_cylinder_DWI.txt"
        N_10_6_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/N_10_6/"+str(conf)+"/dyn_cylinder_DWI.txt"
        paths = [N_10_4_file,N_5_10_4_file, N_10_5_file, N_5_10_5_file, N_10_6_file]
    
        for e, path in enumerate(paths):
            data_dwi = create_data(path)
            x_adc,y_adc,z_adc = get_ADC_value(data_dwi, orientation="separate")
            adcs.append(x_adc)
            orientations.append("x")
            adcs.append(y_adc)
            orientations.append("y")
            adcs.append(z_adc)
            orientations.append("z")
            for i in range(3):
                confs_.append(conf)
                Ns_.append(Ns[e])
    data_adc["ADC"] = adcs
    data_adc["conf"] = confs_
    data_adc["N"] = Ns_
    data_adc["orientation"] = orientations

    fig2 = plt.figure(2)
    ax = sns.lineplot(x="N", y="ADC", hue = "orientation", 
             data=data_adc.loc[data_adc.orientation != "z"] )

    path_ = path.split("/")
    path_ = path_[:-3]
    path_ = "/".join(path_)
    name = path_ + "/ADC_sanity_check.png"

    ax.set_title("ADC wrt number of walkers")
    ax.set_label("number of walkers")
    ax.set_label("ADC [mm²/s]")
    plt.show()
    fig2.savefig(name)

    
def combine_intra_extra_adc(folder_name, rest = True):
    if rest:
        type_ = "rest"
    else :
        type_ = "active"
    extra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(type_)+"/extra/"+str(folder_name)+"/dyn_cylinder_gamma_"+str(type_)+"_extra_DWI.txt"
    intra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(type_)+"/intra/"+str(folder_name)+"/dyn_cylinder_gamma_"+str(type_)+"_intra_DWI.txt"
    bs_extra, extra_dwi = get_total_dwi_array(extra_file)
    bs_intra, intra_dwi = get_total_dwi_array(intra_file)
    fig1 = plt.figure(1)
    plt.plot(bs_extra, extra_dwi, label = "extra")
    plt.plot(bs_intra, intra_dwi, label = "intra")
    dwi_tot = list(map(lambda x,y : x*icvf+(1-icvf)*y, intra_dwi, extra_dwi))
    plt.plot(bs_intra, dwi_tot, label = "total")
    plt.title('DWI Signal')
    plt.ylabel('log(S/So)')
    plt.xlabel('b')
    path_ = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(type_)
    name = path_ + "/DWI_extra_intra_" + str(folder_name)+".png"
    plt.legend()
    fig1.savefig(name)

def combine_active_rest_adc(loc, folders=["N_10_5","N_5_10_5","N_10_6"] , Ns=[1e5, 5e5, 1e6], confs =["conf1", "conf2" , "conf3"] ):

    data_rets_active = pd.DataFrame()
    fs =[] 
    adcs =[] 
    orientations =[] 
    confs_ =[] 
    for e,folder in enumerate(folders): 
        for conf in confs :
            state = "rest"
            path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
            data_dwi = create_data(path)
            x_adc_rest,y_adc_rest,z_adc_rest = get_ADC_value(data_dwi, orientation = "separate")

            state = "active"
            path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
            data_dwi = create_data(path)
            x_adc_act,y_adc_act,z_adc_act = get_ADC_value(data_dwi, orientation = "separate")
            
            adcs.append((x_adc_act-x_adc_rest)*100/x_adc_rest)
            orientations.append("x")
            adcs.append((y_adc_act-y_adc_rest)*100/y_adc_rest)
            orientations.append("y")
            adcs.append((z_adc_act-z_adc_rest)*100/z_adc_rest)
            orientations.append("z")

            for i in range(3):
                fs.append(Ns[e])
                confs_.append(conf)


    data_rets_active["orientation"] = orientations
    data_rets_active["ADC"] = adcs
    data_rets_active["N"] = fs
    data_rets_active["conf"] = confs_
    print(data_rets_active)

    sns.boxplot(data = data_rets_active.loc[data_rets_active.orientation != "z"] , x = "ADC",  y = "orientation", hue = "N")
    plt.show()

state = "rest"
loc = "intra"
#combine_active_rest_adc(loc,folders=["N_10_4"] , Ns=[1e4])
folder = "N_10_4"
conf = "conf1"
path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
            
data = create_data(path)
plot_DWI(data, path)
x,y,z = get_ADC_value(data, orientation = "separate")
mess1 = "Rest: x : "+ str(x) +", y : "+ str(y)
state = "active"
path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
data = create_data(path)
plot_DWI(data, path)
x,y,z = get_ADC_value(data, orientation = "separate")
mess2 = "Active : x : "+ str(x) +", y : "+ str(y)
print(mess1)
print(mess2)
combine_active_rest_adc(loc, folders=["N_10_4"] , Ns=[1e4], confs =["conf1", "conf2", "conf3"] )
