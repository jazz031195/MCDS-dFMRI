# To compare ADC between multiple DWI files with varying N values

import numpy as np
import matplotlib
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

def create_data(dwi_path, b0):
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
    sns.catplot(x="orientations", y="adc [um²/ms]",
             hue="state", col = "loc", kind = "box",
             data=result)

    fig1 = plt.figure(1)
    sns.catplot(x="orientations", y="adc [um²/ms]",
             hue="state", col = "loc", 
             data=result)
    plt.show() 

def get_axons_array(file):
    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 4):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                if (len(spheres)!= 0):
                    axons.append(spheres)
                spheres = []
            

    return axons

def plot_radii(file):
    axons = get_axons_array(file)
    area = 0
    for i in range(len(axons)):
        radius = [] 
        spheres = axons[i]
        for ii in range(len(spheres)):
            radius.append(spheres[ii][3])
        area += np.pi*spheres[ii][3]*spheres[ii][3]
        plt.plot(radius)
    plt.ylabel("[mm]")

    plt.show()


def plot_tortuosity(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
        df = pd.DataFrame()
        df["x[um]"]=xs
        df["y[um]"]=ys
        df["index"] = [i]*len(xs)   
        dfs.append(df)
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="x[um]", y="y[um]")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()

def plot_tortuosity_angle(file):
    axons = get_axons_array(file)
    dfs =[] 
    for j in range(int(len(axons)/100)):
        i = j*50
        xs = [] 
        ys =[]  
        zs =[] 
        spheres = axons[i]
        for ii in range(1,len(spheres)):
            xs.append(spheres[ii][0]-spheres[ii-1][0])
            ys.append(spheres[ii][1]- spheres[ii-1][1])
            zs.append(spheres[ii][2]- spheres[ii-1][2])

        df = pd.DataFrame()
        df["index"] = [i]*len(xs)
        df["theta"]= list(map(lambda x,y:np.arctan2(y,x)/np.pi if np.arctan2(y,x)/np.pi>0 else np.arctan2(y,x)/np.pi+2, xs, ys))
        df["phi"]= list(map(lambda x,y,z:np.arctan(np.sqrt(x**2+y**2)/z)/np.pi if z!= 0 else np.pi/2, xs, ys, zs))  
        dfs.append(df)
        del df
    dfs_ = pd.concat(dfs)
    dfs_ = dfs_.reset_index()
    g = sns.JointGrid(data=dfs_, x="theta", y="phi")
    g.plot_joint(sns.kdeplot,
             fill=True, 
             thresh=0, levels=100, cmap="rocket")
    g.plot_marginals(sns.histplot, color="#03051A", alpha=1, bins=25)
    plt.show()
    

def plot_DWI(intra_rest_path,extra_rest_path, extra_active_path, intra_active_path):
    data_intra_rest = create_data(intra_rest_path)
    data_extra_active = create_data(extra_active_path)
    data_intra_active = create_data(intra_active_path)
    data_extra_rest = create_data(extra_rest_path)
    data_rest = data_intra_rest.add(data_extra_rest)
    data_active = data_intra_active.add(data_extra_active)
    ax = sns.lineplot(data=data_rest, x = "b [ms/um²]", y= "DWI")

    ax = sns.lineplot(data=data_active, x = "b [ms/um²]", y= "DWI")
    plt.show()

def get_total_adc_array(path):
    
    data_dwi = create_data(path)
    data_x = data_dwi.loc[(data_dwi['x'] > 0.0) & (~data_dwi['adc [mm²/s]'].isnull())].sort_values(by = ["b"],ascending=True)
    data_y = data_dwi.loc[(data_dwi['y'] > 0.0) & (~data_dwi['adc [mm²/s]'].isnull())].sort_values(by = ["b"],ascending=True)
    data_z = data_dwi.loc[(data_dwi['z'] > 0.0) & (~data_dwi['adc [mm²/s]'].isnull())].sort_values(by = ["b"],ascending=True)
    bs = list(data_z["b"])
    datas = [data_x,data_y,data_z]

    #x0 = list(dict.fromkeys(list(data_dwi["x"])))[0]
    #data_dwi0 = list(data_dwi.loc[data_dwi.x == x0].sort_values(by = ["b"],ascending=True)["DWI"])
    #x1 = list(dict.fromkeys(list(data_dwi["x"])))[1]
    #data_dwi1 = list(["DWI"])
    #bs = list(data_dwi.loc[data_dwi.x == x1].sort_values(by = ["b"])["b"])
    dwi_tot = list(map(lambda x,y,z : (x+y+z)/3*1000, list(datas[0]["adc [mm²/s]"]),list(datas[1]["adc [mm²/s]"]),list(datas[2]["adc [mm²/s]"])))
    return bs, dwi_tot





def plot_inc_():
    b0 = 0.2
    b1 = 1
    icvf=0.7
    locs = ["intra", "extra"]
    swell = [0.0, 0.3, 0.5, 0.7]
    datas = []
    swell_ = []
    locs_ = []
    trys=[]

    for n in [0,1,2,3,4,5,6,7]:
        for s in swell:
            for l in locs:
                if n == 0:
                    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/icvf_{icvf}_swell_{s}_{l}_DWI.txt"
                else:
                    file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/icvf_{icvf}_swell_{s}_{l}_rep_0{n-1}_DWI.txt"
                print(file)
                if (Path(file)).exists():
                    print("exists")
                    data = get_adc(file, b1, b0)
                    datas.append(data)
                    for i in range(3):
                        swell_.append(s)
                        locs_.append(l)
                        trys.append(n)


    data = pd.concat(datas)
    
    
    data["location"] = locs_
    data["swell_perc"] = swell_
    data["repetition"] = trys
    print(data)
    
    data = data.loc[data["orientations"] != "axial"]
    
    sns.lineplot(data = data.reset_index(), x = "swell_perc",y = "adc [um²/ms]", hue = "location")
    plt.show()
    data = data.groupby(by =["swell_perc","orientations", "location","repetition"] ).mean().reset_index()
    print(data)
    data["weighted_adc"] = list(map(lambda x,y : x*icvf if y == "intra" else x*(1-icvf) , list(data["adc [um²/ms]"] ), list(data["location"])))
    
    data = data.groupby(by =["swell_perc","repetition"] ).mean().reset_index()
    
def combine_intra_extra_adc(folder_name, rest = True):
    # if rest:
    #     type_ = "rest"
    # else :
    #     type_ = "active"


    extra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)+"/extra/_DWI.txt"
    intra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)+"/intra/_rep_00_DWI.txt"
    bs_extra, extra_dwi = get_total_adc_array(extra_file)
    bs_intra, intra_dwi = get_total_adc_array(intra_file)
    bs_intra = [b/1000 for b in bs_intra]
    fig1 = plt.figure(1)
    plt.plot(bs_intra, extra_dwi, label = "extra")
    plt.plot(bs_intra, intra_dwi, label = "intra")
    dwi_tot = list(map(lambda x,y : x*icvf+(1-icvf)*y, intra_dwi, extra_dwi))
    plt.plot(bs_intra, dwi_tot, label = "total")
    plt.title('ADC Signal')
    plt.ylabel('ADC ['r'$\mu$m²/ms]')
    plt.xlabel('b ['r'ms/$\mu$m²]')
    path_ = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)
    name = path_ + "/ADC_extra_intra_.png"
    plt.legend()
    fig1.set_tight_layout(True)
    fig1.savefig(name)

    extra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)+"/extra/_DWI.txt"
    intra_file = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)+"/intra/_rep_00_DWI.txt"
    bs_extra, extra_dwi = get_total_dwi_array(extra_file)
    bs_intra, intra_dwi = get_total_dwi_array(intra_file)
    fig2 = plt.figure(2)
    # extra_dwi = [dwi*1000 for dwi in extra_dwi]
    # intra_dwi = [dwi*1000 for dwi in intra_dwi]
    bs_intra = [b/1000 for b in bs_intra]
    plt.plot(bs_intra, extra_dwi, label = "extra")
    plt.plot(bs_intra, intra_dwi, label = "intra")
    dwi_tot = list(map(lambda x,y : x*icvf+(1-icvf)*y, intra_dwi, extra_dwi)) 
    plt.plot(bs_intra, dwi_tot, label = "total")
    plt.title('DWI Signal')
    plt.ylabel('DWI log(S/S0)')
    plt.xlabel('b ['r'ms/$\mu$m²]')
    path_ = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/"+str(folder_name)
    name = path_ + "/DWI_extra_intra_.png"
    plt.legend()
    fig2.set_tight_layout(True)
    fig2.savefig(name)




def plot_increase():
    b0 = 1
    b1 = 2
    locs = ["intra", "extra"]
    swell = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7]
    datas = []
    swell_ = []
    locs_ = []
    diff_ = []
    w_diff_ = []
    for n in [1, 2]:
        for s in swell:
            for l in locs:
                ref_file = f'/home/localadmin/localdata/juradata/ijelescu/micmap/jasmine/DWI/try{n}/icvf_0.3_swell_0.0_{l}_DWI.txt'
                if (Path(ref_file)).exists():
                    ref_data = get_adc(ref_file, b1, b0)
                    file = f'/home/localadmin/localdata/juradata/ijelescu/micmap/jasmine/DWI/try{n}/icvf_0.3_swell_{s}_{l}_DWI.txt'
                    print(file)
                    if (Path(file)).exists():
                        print("exists")
                        data = get_adc(file, b1, b0)
                        datas.append(get_adc(file, b1, b0))
                        if l == "intra":
                            diff = list((data["adc [um²/ms]"]-ref_data["adc [um²/ms]"])*100/ref_data["adc [um²/ms]"])
                            w_diff_.extend([d*0.3 for d in diff])
                            diff_.extend(diff)
                        else:
                            diff = list((data["adc [um²/ms]"]-ref_data["adc [um²/ms]"])*100/ref_data["adc [um²/ms]"])
                            w_diff_.extend([d*0.7 for d in diff])
                            diff_.extend(diff)
                        for i in range(3):
                            swell_.append(s)
                            locs_.append(l)

    data = pd.concat(datas)
    data["location"] = locs_
    data["swell_perc"] = swell_
    data["adc difference [%]"] = diff_
    data["weighted adc difference [%]"] = w_diff_


    sns.lineplot(data = data.reset_index(), x = "swell_perc",y = "weighted adc difference [%]")
    plt.show()

    sns.lineplot(data = data.reset_index(), hue = "orientations", x = "swell_perc",y = "weighted adc difference [%]")
    plt.show()

    sns.lineplot(data = data.reset_index(), hue = "location", x = "swell_perc",y = "adc difference [%]")
    plt.show()

    data = data.loc[data["swell_perc"] != 0]
    print(data)
    sns.boxplot(data = data.reset_index(),x = "axis",hue = "location",y = "adc difference [%]")
    plt.show()

def plot_const():
    b0 = 0
    b1 = 2
    s = 0.0
    l= "intra"
    datas = []
    repetition = []
    for i in range(5):
        if i == 0:
            file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/icvf_0.3_swell_{s}_{l}_DWI.txt"
        else:
            file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/icvf_0.3_swell_{s}_{l}_rep_0{i-1}_DWI.txt"
        print(file)
        if (Path(file)).exists():
            print("exists")
            data = get_adc(file, b1, b0)
            datas.append(data)
            repetition.extend([i]*len(data))
        
    data = pd.concat(datas)
    data["repetition"] = repetition
    sns.lineplot(data = data.reset_index(), x = "repetition",y = "adc [um²/ms]", hue = "axis")
    plt.show()


    
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


# combine_intra_extra_adc("neurons")
dwi_intra_path = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/intra/_rep_01_DWI.txt"
dwi_intra = create_data(dwi_intra_path)
print(dwi_intra)
dwi_intra["loc"] = ["intra"]*dwi_intra["x"].size
data_x = dwi_intra.loc[(dwi_intra['x'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
data_y = dwi_intra.loc[(dwi_intra['y'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
data_z = dwi_intra.loc[(dwi_intra['z'] > 0.0) ].sort_values(by = ["b [um²/ms]"],ascending=True)
bs = list(data_z["b [um²/ms]"])
datas = [data_x,data_y,data_z]
std = list(map(lambda x,y,z : np.std([x,y,z]), list(datas[0]["log(Sb/So)"]),list(datas[1]["log(Sb/So)"]),list(datas[2]["log(Sb/So)"])))
dwi = dwi_intra

# print(std)
# dwi_extra_path = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/neurons/extra/_DWI.txt"
# dwi_extra = create_data(dwi_extra_path)
# dwi_extra["loc"] = ["extra"]*dwi_extra["x"].size
# data_x = dwi_extra.loc[(dwi_extra['x'] > 0.0) ].sort_values(by = ["b"],ascending=True)
# data_y = dwi_extra.loc[(dwi_extra['y'] > 0.0) ].sort_values(by = ["b"],ascending=True)
# data_z = dwi_extra.loc[(dwi_extra['z'] > 0.0) ].sort_values(by = ["b"],ascending=True)
# bs = list(data_z["b"])
# datas = [data_x,data_y,data_z]
# std = list(map(lambda x,y,z : np.std([x,y,z]), list(datas[0]["log(Sb/So)"]),list(datas[1]["log(Sb/So)"]),list(datas[2]["log(Sb/So)"])))
# print(std)
# dwi = dwi_intra.append(dwi_extra)

import plotly.express as px
fig = px.box(dwi[~dwi['b [um²/ms]'].isnull()], x='b [um²/ms]', y='Sb/So', color='loc')
fig.show()
# intra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/intra_active__DWI.txt"
# extra_active_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/extra_active__DWI.txt"
# intra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/intra__DWI.txt"
# extra_rest_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/axons/extra__DWI.txt"               
# data2 = assemble_data(intra_active_path, intra_rest_path, extra_active_path, extra_rest_path)

# print(data)

# plot(data2)