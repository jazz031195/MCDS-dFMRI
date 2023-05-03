import numpy as np
import os
import random
import warnings
warnings.filterwarnings("ignore")

cur_path = os.getcwd()

# Creates copies of the axon list while modifying the number of active neurons. 
# Some axons that could previously swell will become static

def get_axons_array(file):
    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 4):
                # add the activity
                spheres.append(float(line.split(' ')[-1]))
            else :
                if (len(spheres)!= 0):
                    axons.append(np.mean(spheres))
                spheres = []
    print("Percentage swelling in "+str(file)+ " is :" + str(np.sum(axons)/len(axons)))

    return axons

def make_new(perc_swell, axons_file):

    axons = get_axons_array(axons_file)
    if (perc_swell> np.sum(axons)/len(axons)):
        print("Select lower perc_swell")
        return
    # difference between active axons initially and after
    diff = int(np.sum(axons)-int(len(axons)*perc_swell))
    # create list that shows which active axons are active in new file
    lst = [1]*int(np.sum(axons))
    lst[:diff] = [0]*diff
    random.shuffle(lst)
    axon_activity = axons
    n = 0
    for e,i in enumerate(axon_activity) :
        if i == 1:
            axon_activity[e] = lst[n]
            n += 1

    # create new file
    path_to_new_file = ("/").join(axons_file.split("/")[:-1])
    new_file = axons_file.split("/")[-1].split("_")
    new_file = ("_").join(new_file[:3])+"_"+str(perc_swell)+"_"+("_").join(new_file[4:])
    new_file = path_to_new_file +"/"+ new_file

    with open(axons_file,'r') as firstfile, open(new_file,'w') as secondfile:
        axon_nbr = 0

        for e,line in enumerate(firstfile):
                if (len(line.split(' ')) > 4):
                    activity = axon_activity[axon_nbr]
                    new_line = line.split(" ")
                    new_line[-1] = str(int(activity))
                    new_line = (" ").join(new_line)
                    secondfile.write(new_line)
                    secondfile.write("\n")
                    
                else:
                    if e > 5:
                        axon_nbr += 1
                    # write content to second file
                    secondfile.write(line)
    return new_file


for i in [0.0]:
    perc = 0.3
    axons_file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/axons/Substrates/icvf_0.3_swell_{perc}_gamma_distributed_axon_list.txt"
    new_file = make_new(i, axons_file)
    get_axons_array(new_file)