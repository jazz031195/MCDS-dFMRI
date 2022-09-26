import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import os

file_name = "cylinder_gamma_packing_test_0.traj.txt"
path = '/MCDC_Simulator_public-master/MCDC_Simulator_public-master/instructions/demos/output/' + str(file_name)

cur_path = os.path.dirname(__file__)

new_path = os.path.relpath(path, cur_path)
print(new_path)
file = open(new_path, mode="r")
#contents = file.read()