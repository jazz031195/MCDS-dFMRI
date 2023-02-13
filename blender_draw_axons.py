import bpy

def get_spheres_array(file):

    axons = []
    spheres = []

    with open(file) as f:
        for line in f.readlines():
            if (len(line.split(' ')) > 3):
                spheres.append([float(i) for i in line.split(' ')[:]])
            else :
                axons.append(spheres)
                spheres = []

    return axons

file = "/home/localadmin/Documents/MCDS_code/MCDS-dFMRI/MCDC_Simulator_public-master/instructions/demos/output/axons/_gamma_distributed_axon_list_.txt"
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

axons = get_spheres_array(file)
for i in range(len(axons)):
    spheres = axons[i]
    for ii in range(len(spheres)):
        r = spheres[ii][3]
        x = spheres[ii][0]
        y = spheres[ii][1]
        z = spheres[ii][2]

        bpy.ops.mesh.primitive_uv_sphere_add(radius=r, enter_editmode=False, align='WORLD', location=(x, y, z), scale=(100, 100, 100))



