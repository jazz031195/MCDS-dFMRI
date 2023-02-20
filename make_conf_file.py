
import os
import sys
import getopt

cur_path = os.getcwd()

def write_conf_file_create_axons(N, T, exp_prefix, icvf , num_axons, c2):
    """
    Overwrites the configuration file gammaDistributedAxons.conf with the desired parameters
    Should be used in case we want to grow axons in a new substrate
    N : number of water particles
    T : number of time steps in simulation
    exp_prefix : name of file
    icvf : Intracompartement volume fraction
    num_axons : number of axons
    c2 : measure of ODF 
    location : can be intra or extra

    """
    
    lines =[]; 
    lines.append("N "+ str(N))
    lines.append("T "+ str(T))
    lines.append("duration 0.067")
    lines.append("diffusivity 2.5e-9")
    lines.append("active_state false")
    lines.append("exp_prefix /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/" +exp_prefix)
    lines.append("scheme_file /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme")
    lines.append("scale_from_stu 1")
    lines.append("write_txt 1")
    lines.append("write_bin 0")
    lines.append("write_traj_file 0")
    lines.append("<obstacle>")
    lines.append("<axon_gamma_packing>")
    lines.append("alpha 5")
    lines.append("beta 0.1")
    lines.append("icvf "+str(icvf))
    lines.append("num_axons "+str(num_axons))
    lines.append("percentage_dynamic_axons 0.3")
    lines.append("percentage_increase_of_volume 0.01")
    lines.append("c2 "+str(c2))
    lines.append("tortuous true")
    lines.append("</axon_gamma_packing>")
    lines.append("</obstacle>")
    lines.append("ini_walkers_pos extra")
    lines.append("num_process 10")
    lines.append("<END>")

    name = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

    f = open(name, "w")
    for l in lines :
        f.write(str(l))
        f.write("\n")
    f.close()

def write_conf_file(N, T, activation, exp_prefix, location, axons_list_prefix):

    """
    Overwrites the configuration file gammaDistributedAxons.conf with the desired parameters
    Should be used in case we want to use a substrate from an existing file
    N : number of water particles
    T : number of time steps in simulation
    activation : can be active (axons swell) or rest
    exp_prefix : name of file
    location : can be intra or extra
    axons_list_prefix : name of file with list of axons

    """
    lines =[]; 
    lines.append("N "+ str(N))
    lines.append("T "+ str(T))
    lines.append("duration 0.067")
    lines.append("diffusivity 2.5e-9")
    lines.append("active_state " + str(activation))
    lines.append("exp_prefix /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/" +exp_prefix)
    lines.append("scheme_file /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/scheme_files/PGSE_sample_scheme_new.scheme")
    lines.append("scale_from_stu 1")
    lines.append("write_txt 1")
    lines.append("write_bin 0")
    lines.append("write_traj_file 0")
    lines.append("<obstacle>")
    lines.append("axon_list /home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/"+str(axons_list_prefix)+"_gamma_distributed_axon_list.txt")
    lines.append("</obstacle>")
    lines.append("ini_walkers_pos "+ str(location))
    lines.append("num_process 10")
    lines.append("<END>")

    name = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/conf_file_examples/gammaDistributedAxons.conf"

    f = open(name, "w")
    for l in lines :
        f.write(str(l))
        f.write("\n")
    f.close()



def get_variables():

    """
    Obtains variables through input in bash script

    """

    nbr_axons  = None
    T = 1000
    loc = None 
    state = None
    icvf = None
    concentration = 1000
    create_substrate_ = "true"

    argv = sys.argv[1:]
  
    try:
        opts, args = getopt.getopt(argv, ":a:t:l:s:i:c:x:")
      
    except:
        print("Error")
  
    for opt, arg in opts:
        if opt in ['-a']:
            # number of axons 
            nbr_axons = arg
        elif opt in ['-t']:
            # time   
            T = arg
        elif opt in ['-l']:
            # location : extra or intra 
            loc = arg
        elif opt in ['-s']:
            # state : active or rest
            state = arg
        elif opt in ['-i']:
            # icvf 
            icvf = arg
        elif opt in ['-c']:
            # c2 
            c2 = float(arg)
        elif opt in ['-x']:
            create_substrate_ = arg
        
    if (create_substrate_ == "true"):
        create_substrate = True
    elif (create_substrate_ == "false"):
        create_substrate = False

    # number of particles is set as a function on the compartment and icvf choices 
    if (loc == "intra") :
        N = int(float(icvf)*concentration)
    else :
        N = int((1-float(icvf))*concentration)

        
        
    return N, nbr_axons, T, loc, state,icvf, create_substrate, c2
        
      

N, nbr_axons, T, location, activation, icvf, create_substrate, c2 = get_variables()
axons_list_prefix = "nbr_axons_"+str(nbr_axons)+"_icvf_"+str(icvf)
if(create_substrate):
    write_conf_file_create_axons(1, 1, axons_list_prefix, icvf , nbr_axons, c2)
else:
    exp_prefix = "nbr_axons_"+str(nbr_axons)+"_"+str(location)+"_"+str(activation)+"_icvf_"+str(icvf)
    write_conf_file(N, T, activation, exp_prefix, location, axons_list_prefix)
