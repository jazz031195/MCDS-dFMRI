# To compare ADC between multiple DWI files with varying N values

from pathlib import Path
import sys
import getopt

mcdc_dir = Path() / "MCDC_Simulator_public-master"
dynamic_cylinders = mcdc_dir / "instructions/demos/output/dynamic_cylinders"
conf_file_examples = mcdc_dir / "docs/conf_file_examples"
scheme_files = mcdc_dir / "docs/scheme_files"

def get_variables():
    N = None 
    conf = None 
    loc = None 
    state = None
    folder = None
    create_substrate = 0

    argv = sys.argv[1:]
  
    try:
        opts, args = getopt.getopt(argv, ":n:c:l:s:f:d:")
      
    except:
        print("Error")
  
    for opt, arg in opts:
        if opt in ['-n']:
            N = arg
        elif opt in ['-c']:
            conf = arg
        elif opt in ['-l']:
            loc = arg
        elif opt in ['-s']:
            state = arg
        elif opt in ['-f']:
            folder = arg
        elif opt in ['-d']:
            create_substrate = int(arg)
        
        if create_substrate == 1:
            create_substrate = True
        else :
            create_substrate = False

        
        
    return N, conf, loc, folder, state, create_substrate
        
      

def create_lines(N, conf, loc, folder, state, create_substrate):
    if create_substrate:
        conf_file = conf_file_examples / "create_substrate.conf"

    else:
        conf_file = conf_file_examples / "read_from_substrate.conf"


    if state == "rest":
        active_state = "false"
    elif state == "active":
        active_state = "true"

    new_lines =[] 
    with open(conf_file) as f:
        for line in f.readlines():
            e = line.split(' ')
            if e[0] == "N":
                l = f"N {N}"
            elif e[0] == "active_state":
                l = f"active_state {active_state}"
            elif e[0] == "exp_prefix":
                if not create_substrate:
                    prefix = dynamic_cylinders / state / loc / folder / conf
                else:
                    prefix = dynamic_cylinders
                if not prefix.exists():
                    prefix.mkdir()
                l = f"exp_prefix {prefix / 'dyn_cylinder'}"

            # Make sure to have a compatible path also for scheme_file
            elif e[0] == "scheme_file":
                l = f"scheme_file {scheme_files / Path(e[1]).name}"
           
            elif e[0] == "ini_walkers_pos":
                l = f"ini_walkers_pos {loc}"
            
            elif e[0] == "dyn_cylinders_list":
                if not create_substrate:
                    l = f"dyn_cylinder_list {dynamic_cylinders / conf}.txt"
                else :
                    l = line
            else :
                l = line

            new_lines.append(l)

    return new_lines

def write_conf_file(lines):

    name = conf_file_examples / "gammaDistributedCylinders.conf"

    f = open(name, "w")
    for l in lines :
        f.write(str(l))
        f.write("\n")
    f.close()

N, conf, loc, folder, state, create_substrate = get_variables()
lines = create_lines(N, conf, loc, folder, state, create_substrate)
write_conf_file(lines)
