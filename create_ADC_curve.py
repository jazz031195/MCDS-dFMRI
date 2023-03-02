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
        b0 = list(datas[i]["b"])[0]
        Sb0 = list(datas[i].loc[datas[i]["b"]== b0]["DWI"])[0]
        signal = list(map(lambda Sb : np.log(Sb/Sb0), list(datas[i]["DWI"])))
        adc = list(map(lambda b,Sb : -np.log(Sb/Sb0)/(b-b0) if b!= b0 else np.nan, list(datas[i]["b"]),list(datas[i]["DWI"])))
        datas[i]["log(Sb/So)"] = signal
        datas[i]["adc [mm²/s]"] = adc

    data_dwi = pd.concat(datas)

    return data_dwi

def get_adc(dwi_path):
    data = create_data(dwi_path)
    axes = ["x","y","z"]
    orientations = ["radial","radial","axial"]
    adcs = []
    for a in axes:
        ax_data = data.loc[data[a]>0]
        b1 = list(ax_data["b"])[-1]
        ax_data = ax_data.loc[ax_data["b"] == b1]
        adcs.append(list(ax_data["adc [mm²/s]"])[0])
    new_data = pd.DataFrame(columns = ["axis", "adc [mm²/s]"])
    new_data["orientations"] = orientations
    new_data["axis"] = axes
    new_data["adc [mm²/s]"] = adcs
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
    sns.catplot(x="orientations", y="adc [mm²/s]",
             hue="state", col = "loc", kind = "box",
             data=result)

    fig1 = plt.figure(1)
    sns.catplot(x="orientations", y="adc [mm²/s]",
             hue="state", col = "loc", 
             data=result)
    plt.show() 



def plot_DWI(data, path):
    data["orientation"] =  list(data["x"]).copy()
    data["orientation"] = list(map(lambda x,y : "x" if x == 1 else y,list(data["x"]), list(data["orientation"])))
    data["orientation"] = list(map(lambda x,y : "y" if x == 1 else y,list(data["y"]), list(data["orientation"])))
    data["orientation"] = list(map(lambda x,y : "z" if x == 1 else y,list(data["z"]), list(data["orientation"])))
    fig1 = plt.figure(1)
    sns.lineplot(x="b", y="log(Sb/So)",
             hue="orientation",
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
        return np.mean([adcs[0],adcs[1],2.5e-3])
        #return ((adcs[0]+adcs[1])/2)
    if orientation == "radial":
        return ((adcs[0]+adcs[1])/2)
    if orientation == "z":
        return (adcs[2])
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
            
            if Path(path).exists():

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
             data=data_adc.loc[data_adc.orientation == "z"] )
    

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



def combine_active_rest_adc(folders=["N_10_5","N_5_10_5","N_10_6"] , Ns=[1e5, 5e5, 1e6], confs =["conf1", "conf2" , "conf3"] ):

    data_rets_active = pd.DataFrame()
    fs =[]
    adcs =[]
    orientations =[] 
    confs_ =[] 
    locs_ = []
    locs = ["extra", "intra"]
    for loc in locs:
        for e,folder in enumerate(folders): 
            for conf in confs :
                state = "rest"
                path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
                if Path(path).exists():
                    data_dwi = create_data(path)
                    x_adc_rest,y_adc_rest,z_adc_rest = get_ADC_value(data_dwi, orientation = "separate")

                    state = "active"
                    path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI.txt"
                    if Path(path).exists():
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
                            locs_.append(loc)

    data_rets_active["orientation"] = orientations
    data_rets_active["% ADC change from rest to active"] = adcs
    data_rets_active["N"] = fs
    data_rets_active["conf"] = confs_
    data_rets_active["location"] = locs_
    print(data_rets_active)

    sns.stripplot(data = data_rets_active , x = "% ADC change from rest to active",  y = "orientation", hue = "conf")
    plt.show()

def final_plot_adc():
    folder = "N_10_6"
    conf = "conf2"
    states = ["active", "rest"]
    locs =["intra", "extra"] 
    df = pd.DataFrame()
    adcs =[]
    locs_ = []
    states_ = []
    orientations =[] 
    trials_ = []
    for loc in locs :
        for state in states:
            for i in range(8):
                path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
                data = create_data(path)
                x_adc_act,y_adc_act,z_adc_act = get_ADC_value(data, orientation = "separate")
                adcs.append((x_adc_act+y_adc_act)/2)
                orientations.append("radial")
                adcs.append(z_adc_act)
                orientations.append("z")

                for ii in range(2):
                    locs_.append(loc)
                    states_.append(state)
                    trials_.append(i)
    df["orientation"] = orientations
    df["ADC [mm²/s]"] = adcs
    df["location"] = locs_
    df["state"] = states_
    df = df.loc[df.orientation == "radial"]
    print(df)
    plt.figure(figsize=(10,8))
    #g = sns.FacetGrid(df.loc[df.orientation == "radial"], y= "ADC [mm²/s]", row = "orientation", hue="state",sharex=False, hue_order= [ "rest", "active"])
    #g.map(sns.catplot, "ADC [mm²/s]", x= "location", alpha=.7)
    #g.set_xticklabels(fontsize=5)
    #g.add_legend()
    g = sns.catplot(data = df, y= "ADC [mm²/s]", x= "state",  hue = "state", col = "location", hue_order=["rest","active"], sharey=False)
    #plt.show()
    for axs in g.axes:
        for ax in axs:
            ax.set_yticklabels(ax.get_yticklabels(),fontsize= 12)
            ax.set_xticklabels(ax.get_xticklabels(),fontsize= 12)

    plt.savefig("adc_compare.png")

def stats():
    folder = "N_10_6"
    conf = "conf2"
    states = ["active", "rest"]
    df = pd.DataFrame()
    rs =[]
    zs = []
    states_ = []

    for state in states :
        for i in range(5):
            loc = "intra"
            path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
            data = create_data(path)
            x_adc_intra,y_adc_intra,z_adc_intra = get_ADC_value(data, orientation = "separate")

            loc = "extra"
            path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
            data = create_data(path)
            x_adc_extra,y_adc_extra,z_adc_extra = get_ADC_value(data, orientation = "separate")
            icvf = 0.7
            if state == "active":
                icvf *= (1+0.01*0.3)
            x_adc = (1-icvf)*x_adc_extra+ icvf*x_adc_intra
            y_adc = (1-icvf)*y_adc_extra+ icvf*y_adc_intra
            z_adc = (1-icvf)*z_adc_extra+ icvf*z_adc_intra

            rs.append((x_adc+y_adc)/2)
            zs.append(z_adc)

            states_.append(state)
            
    df["rADC"] = rs
    df["zADC"] = zs
    df["state"] = [1 if s == "active" else 0 for s in states_]

    
    y = df.state
    X = df[["rADC","zADC"]]

    X2 = sm.add_constant(X)
    est = sm.OLS(y, X2)
    est2 = est.fit()
    print(est2.summary())





def get_adc_diff_rest_active(folder, conf):

    adcs = []
    df_adc = pd.DataFrame()
    states = []
    trial = []
    adcs_r = []
    rests = []
    directions = []

    plot_types = ["all","radial"]
    
    for i in range(8):
        icvf = 0.7
        state = "rest"
        loc = "intra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_intra = get_ADC_value(data, orientation = plot_types[0])
        loc = "extra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_extra = get_ADC_value(data, orientation = plot_types[0])
        adc_rest = (1-icvf)*adc_extra + icvf * adc_intra
        adcs.append(adc_rest)
        rests.append(adc_rest)
        directions.append("all directions")

        loc = "intra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_intra = get_ADC_value(data, orientation = plot_types[1])
        loc = "extra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_extra = get_ADC_value(data, orientation = plot_types[1])
        adc_rest = (1-icvf)*adc_extra + icvf * adc_intra
        adcs.append(adc_rest)
        directions.append("radial")

        for ii in range(2):
            states.append("rest")
            trial.append(i)

        icvf *= (1+ 0.3*0.01) 
        state = "active"
        loc = "intra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_intra = get_ADC_value(data, orientation = plot_types[0])
        loc = "extra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_extra = get_ADC_value(data, orientation = plot_types[0])
        adc_active= (1-icvf)*adc_extra + icvf * adc_intra
        adcs.append(adc_active)
        rests.append(adc_active)
        directions.append("all directions")

        loc = "intra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_intra = get_ADC_value(data, orientation = plot_types[1])
        loc = "extra"
        path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/dynamic_cylinders/"+str(state)+"/"+str(loc)+"/"+str(folder)+"/"+str(conf)+"/dyn_cylinder_DWI_"+str(i)+".txt"
        data = create_data(path)
        adc_extra = get_ADC_value(data, orientation = plot_types[1])
        adc_active= (1-icvf)*adc_extra + icvf * adc_intra
        adcs.append(adc_active)
        directions.append("radial")

        for ii in range(2):
            states.append("task")
            trial.append(i)


    df_adc["state"] = states
    df_adc["ADC [mm²/s]"] = adcs
    df_adc["direction"] = directions
    df_adc["Relative adc"] = (df_adc["ADC [mm²/s]"]/np.mean(rests))
    df_adc["trial"] = trial


    g = sns.catplot(data = df_adc, x = "state", y = "ADC [mm²/s]", col = "direction", sharey=False)

    for axs in g.axes:
        for ax in axs:
            ax.set_yticklabels(ax.get_yticklabels(),fontsize= 15)
            ax.set_xticklabels(ax.get_xticklabels(),fontsize= 15)


    plt.savefig("ADC_tot.png")

    df = df_adc.groupby(by = ["state"]).mean()
    active = list(df["Relative adc"])[0]
    rest = list(df["Relative adc"])[1]
    print("% Change in ADC from rest to active : ", (active-rest)*100/rest)
    print("Change in ADC from rest to active : ", (active-rest))


#dwi_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/DWI_intra_rest.txt"

#data = create_data(dwi_path)

intra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/nbr_axons_50_intra_active_icvf_0.5_DWI.txt"
extra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/nbr_axons_50_extra_active_icvf_0.5_DWI.txt"
intra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/nbr_axons_50_intra_rest_icvf_0.5_DWI.txt"
extra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/nbr_axons_50_extra_rest_icvf_0.5_DWI.txt"               
data2 = assemble_data(intra_active_path, intra_rest_path, extra_active_path, extra_rest_path)

#print(data)

plot(data2)