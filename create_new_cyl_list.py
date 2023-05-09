import numpy as np
import os
import random
import warnings
warnings.filterwarnings("ignore")

cur_path = os.getcwd()

# Creates copies of the axon list while modifying the number of active neurons. 
# Some axons that could previously swell will become static

def get_cyls_array(file):

    cyls = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 4):
                # add the activity
                cyls.append(float(line.split(' ')[-1]))

    swelling_cyls = np.sum(cyls)
    print("Percentage swelling in "+str(file)+ " is :" + str(swelling_cyls/len(cyls)))

    return cyls

def make_new(old_perc, new_perc, cyls_file):

    cyls = get_cyls_array(cyls_file)
    if (new_perc> old_perc):
        print("Select lower perc_swell")
        return
    # number of swelling cylinders
    swell_cyl = np.sum(cyls)
    # difference between active axons initially and after
    diff = int(swell_cyl-new_perc*len(cyls))
    # create list that shows which active axons are active in file
    lst = [1]*int(swell_cyl)
    # reset to 0 some of the active cylinders 
    lst[:diff] = [0]*diff
    # shuffle list
    random.shuffle(lst)
    axon_activity = cyls
    # assign new activity whenever activity = 1
    n = 0
    for e,i in enumerate(axon_activity) :
        if i > 0:
            axon_activity[e] = lst[n]
            n += 1

    # create new file
    path_to_new_file = ("/").join(cyls_file.split("/")[:-1])
    new_file = cyls_file.split("/")[-1].split("_")
    new_file = ("_").join(new_file[:3])+"_"+str(new_perc)+"_"+("_").join(new_file[4:])
    new_file = path_to_new_file +"/"+ new_file

    with open(axons_file,'r') as firstfile, open(new_file,'w') as secondfile:
        cyl_nbr = 0

        for e,line in enumerate(firstfile):
                if (len(line.split(' ')) > 4):
                    activity = axon_activity[cyl_nbr]
                    new_line = line.split(" ")
                    new_line[-1] = str(int(activity))
                    new_line = (" ").join(new_line)
                    secondfile.write(new_line)
                    secondfile.write("\n")
                    cyl_nbr += 1
                    
                else:
                    # write content to second file
                    secondfile.write(line)
    return new_file

lst = [0.7, 0.6, 0.5, 0.4, 0.3, 0.0]
for i,v in enumerate(lst):
    if i > 0:
        last_perc = lst[i-1]
    else:
        last_perc = 1.0
    axons_file = cur_path + f"/MCDC_Simulator_public-master/instructions/demos/output/cylinders/Substrates/icvf_0.7_swell_{last_perc}_gamma_distributed_dyn_cylinder_list.txt"
    new_file = make_new(last_perc, v,  axons_file)
