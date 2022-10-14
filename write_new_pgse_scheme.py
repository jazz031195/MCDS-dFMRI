import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

cur_path = os.getcwd()
scheme_file_path = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme.scheme"
#delta = 10 # ms
#Delta = 50 #ms
b_values = [0, 1] # ms/um²
giro = 267.51525e3 #Gyromagnetic radio given in rad/(ms*T)
filename_0 = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_b0.scheme"
filename_1 = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_b1.scheme"

def read_scheme(b_value):
    with open(scheme_file_path) as f:
        all_lines = []
        for line in f.readlines():
            if len(line.split(' ')) > 3:
                new_line = []
                Delta = float(line.split(' ')[4])*1000
                delta = float(line.split(' ')[5])*1000
                for e, element in enumerate(line.split(' ')):
                    if e ==3 :
                        G = np.sqrt(b_value/(Delta-(delta/3)))/(giro*delta)
                        new_line.append(str(round(float(G)*1e6,3)))# um to m
                    #elif e == 4 :
                    #    new_line.append(str(float(Delta)*0.001))# ms to s
                    #elif e == 5:
                    #    new_line.append(str(float(delta)*0.001))# ms to s
                    else: 
                        new_line.append(element)
                all_lines.append(new_line)
            
    return all_lines

def write_scheme(all_lines, b_value):
    if b_value == 0:
        name = filename_0
    elif b_value ==1 :
        name = filename_1
    f = open(name, "a")
    for i,line in enumerate(all_lines):
        if i == 0:
            f.write("VERSION: JASMINE")
            f.write("\n")
        f.write(" ".join(line))
    f.close()

# main
for b in b_values:
    all_lines = read_scheme(b)
    write_scheme(all_lines, b)