
import os
import sys
import getopt

cur_path = os.getcwd()

path = "/scratch/jnguyend/Simulation/MCDC_Simulator_public-master"

def write_conf_file_create_axons(N, T, exp_prefix, icvf , num_axons, c2, i, swell):
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
    lines.append(f"exp_prefix {path}/instructions/demos/output/axons/" +exp_prefix)
    lines.append(f"scheme_file {path}/docs/scheme_files/PGSE_sample_scheme_new.scheme")
    lines.append("scale_from_stu 1")
    lines.append("write_txt 1")
    lines.append("write_bin 0")
    lines.append("write_traj_file 0")
    lines.append("<obstacle>")
    lines.append("<axon_gamma_packing>")
    lines.append("alpha 0.0")
    lines.append("beta 0.0")
    lines.append("icvf "+str(icvf))
    lines.append("num_axons "+str(num_axons))
    lines.append("percentage_dynamic_axons "+str(swell))
    lines.append("percentage_increase_of_volume 0.01")
    lines.append("c2 "+str(c2))
    lines.append("tortuous true")
    lines.append("min_radius 0.1")
    lines.append("</axon_gamma_packing>")
    lines.append("</obstacle>")
    lines.append("<voxels>")
    lines.append("0 0 0")
    lines.append("0 0 0")
    lines.append("</voxels>")
    lines.append("ini_walkers_pos extra")
    lines.append("num_process 10")
    lines.append("<END>")

    name = f"{path}/docs/conf_file_examples/gammaDistributedAxons.conf"

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
    lines.append(f"exp_prefix {path}/instructions/demos/output/axons/" +exp_prefix)
    lines.append(f"scheme_file {path}/docs/scheme_files/PGSE_sample_scheme_new.scheme")
    lines.append("scale_from_stu 1")
    lines.append("write_txt 1")
    lines.append("write_bin 0")
    lines.append("write_traj_file 0")
    lines.append("<obstacle>")
    lines.append(f"axon_list {path}/instructions/demos/output/axons/Substrates/"+str(axons_list_prefix)+".txt")
    lines.append("</obstacle>")
    lines.append("ini_walkers_pos "+ str(location))
    lines.append("num_process 20")
    lines.append("<END>")

    name = f"{path}/docs/conf_file_examples/gammaDistributedAxons.conf"

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
    concentration = 100000
    create_substrate_ = "true"
    g = None
    swell = None

    argv = sys.argv[1:]
  
    try:
        opts, args = getopt.getopt(argv, ":a:t:l:s:i:c:x:g:S:")
      
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
        elif opt in ['-g']:
            g = arg
        elif opt in ['-S']:
            # percentage of axons swelling
            swell = arg
        
    if (create_substrate_ == "true"):
        create_substrate = True
    elif (create_substrate_ == "false"):
        create_substrate = False

    # number of particles is set as a function on the compartment and icvf choices 
    if (loc == "intra") :
        N = int(float(icvf)*concentration)
    else :
        N = int((1-float(icvf))*concentration)

        
        
    return N, nbr_axons, T, loc, state,icvf, create_substrate, c2, g, swell
        
      

N, nbr_axons, T, location, activation, icvf, create_substrate, c2, g, swell = get_variables()
#axons_list_prefix = "nbr_axons_"+str(nbr_axons)+"_icvf_"+str(icvf)
axons_list_prefix = f'icvf_{icvf}_swell_{swell}_gamma_distributed_axon_list'
print(axons_list_prefix )
if(create_substrate):
    write_conf_file_create_axons(1, 1, axons_list_prefix, icvf , nbr_axons, c2, g, swell)
else:
    exp_prefix = f'icvf_{icvf}_swell_{swell}_{location}'
    write_conf_file(N, T, activation, exp_prefix, location, axons_list_prefix)
