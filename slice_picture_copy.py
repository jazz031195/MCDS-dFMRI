import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os
from scipy.linalg import norm

cur_path = os.getcwd()
cyl_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_10_5/dyn_cylinder_gamma_rest_extra_gamma_distributed_dyn_cylinder_list.txt"
new_cyl_file_path = cur_path + "/MCDC_Simulator_public-master/instructions/demos/output/sanity_check/N_10_5/dyn_cylinder_gamma_rest_extra_gamma_distributed_dyn_cylinder_list_new.txt"

with open(cyl_file_path) as f:
    list_ = []
    for line in f.readlines():
        if len(line.split(' ')) > 3:
            l = (line[:-2]).split(' ')
            l.append(line[-2])
            print(l)
        else :
            l = []
            l.append(line)
        list_.append(l)
with open(new_cyl_file_path, "w") as file:
    for l in list_:
        if len(l) >3 :
            file.write(" ".join(l))
            file.write("\n")
        else :
            file.write(l[0])
        