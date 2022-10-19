import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

cur_path = os.getcwd()
delta = 16.5e-3 # ms
Delta = 50e-3 #ms
b_values = [0, 200, 500, 1000] # ms/um²
TE = 0.1
giro = 2.6751525e8 #Gyromagnetic radio given in rad/(ms*T)
filename = cur_path + "/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme"

def read_scheme(b_value):
    with open(filename) as f:
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

def write_scheme(vectors, b_values, delta, Delta, TE):
    name = filename
    f = open(name, "a")
    f.write("VERSION: STEJSKALTANNER")
    f.write("\n")
    for vector in vectors :
        for b in b_values :
            G = round(np.sqrt(b/(pow(giro*delta,2)*(Delta-delta/3)))*1e3,3)
            for v in vector:
                f.write(str(v))
                f.write(" ")
            f.write(str(G))
            f.write(" ")
            f.write(str(Delta))
            f.write(" ")
            f.write(str(delta))
            f.write(" ")
            f.write(str(TE))
            f.write("\n")
    f.close()

# main
v1 = np.random.rand(2)
v2 = np.random.randn(2)
v1 /= np.linalg.norm(v1) 
v2 -= v2.dot(v1) * v1
v2 /= np.linalg.norm(v2) 

v1 = [round(v,3) for v in v1]
v2 = [round(v,3) for v in v2]
v1.append(0.0)
v2.append(0.0)
vectors = [v1,v2]
write_scheme(vectors, b_values, delta, Delta, TE)