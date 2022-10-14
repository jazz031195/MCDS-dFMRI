import numpy as np
import matplotlib.pyplot as plt
import os
cur_path = os.getcwd()
scheme_file_path_b0 = cur_path + "/MCDC_Simulator_public-master/docs/instructions/demos/output/dyn_cylinder_gamma_packing_test_b0_DWI.txt"
scheme_file_path_b1 = cur_path + "/MCDC_Simulator_public-master/docs/instructions/demos/output/dyn_cylinder_gamma_packing_test_b1_DWI.txt"

def read_dwi_file(b_value):
    if b_value == 0:
        file = scheme_file_path_b0
    elif b_value == 1:
        file = scheme_file_path_b1
    with open(file) as f:
        dwi_signal = []
        for line in file.readlines():
            dwi_signal.append(float(line))
    return np.array(dwi_signal)

def create_ADC(dwi_signal0, dwi_signal1, b_value):
    dv = np.divide(dwi_signal1, dwi_signal0)
    adc = np.array(list(map(lambda x: -np.log(x)/b_value, dv)))
    plt.plot(range(adc),adc)
    title_ = "Apparent diffusion coefficient, b =" + str(b_value)
    plt.title(title_)
    plt.show()
    return adc

dwi_signal0 = read_dwi_file(0)
dwi_signal1 = read_dwi_file(1)
create_ADC(dwi_signal0, dwi_signal1, 1)