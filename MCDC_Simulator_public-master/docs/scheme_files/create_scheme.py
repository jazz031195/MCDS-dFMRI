import numpy as np

# Specify the path to the file
file_path = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/scheme_files/samples.txt"

# Initialize empty lists to store the gradient components
dir = []

# Read the file line by line and extract the gradient components
with open(file_path, "r") as file:
    lines = file.readlines()
    for line in lines:
        # Skip comment lines and empty lines
        if line.strip() == "" or line.startswith("#"):
            continue
        # Split the line into columns
        columns = line.strip().split("\t")
        # Extract the ux, uy, and uz values from the columns
        ux = float(columns[1])
        uy = float(columns[2])
        uz = float(columns[3])
        # Append the values to the respective lists
        dir.append([ux, uy, uz])
        print(np.linalg.norm([ux, uy, uz]))

dir = np.array(dir)
# Now you have the extracted gradient components in the lists ux_values, uy_values, and uz_values
print("dir:", dir)

# Specify the path to the file
file_path = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/docs/scheme_files/PGSE_21_dir.scheme"

G = [0.0, 0.015, 0.034, 0.048, 0.059, 0.107]
Delta = 0.05 #s
delta = 0.0165 #s
TE    = 0.067 #s
# Open the file in write mode
# If the file doesn't exist, it will be created. If it exists, its contents will be overwritten.
with open(file_path, "w") as file:
    file.write("VERSION: STEJSKALTANNER \n")
    for d in dir:
        for g in G:
            # Write the data to the file
            file.write(f"{d[0]} {d[1]} {d[2]} {g} {Delta} {delta} {TE}\n")