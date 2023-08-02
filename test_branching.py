import numpy as np
import matplotlib.pyplot as plt

def rotateDirection(direction, angle):
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    v1 = direction[:2].dot(rotation_matrix)
    v1 = np.append(v1, 0)
    v1 = v1 / np.linalg.norm(v1)

    return v1


parent_dir       = np.array([1, 0, 0])
bifurcationAngle = (np.pi/4 - np.pi/16) / 2
children_dir = []
for i in range(2):
    rotatedDir = rotateDirection(parent_dir, bifurcationAngle)
    children_dir.append(rotatedDir)
    bifurcationAngle = -bifurcationAngle

list_spheres = np.zeros((1, 3))
center = np.zeros(3)
radius = 0.5e-3
for i in range(1, 21):
    list_spheres = np.vstack((list_spheres, center + parent_dir * radius / 4 * i))

last_sphere = list_spheres[-1]
for j in range(2):
    for i in range(1, 21):
        list_spheres = np.vstack((list_spheres, last_sphere + (i+1) * children_dir[j] * radius / 4))
fig, axs = plt.subplots(1, 1)
for i in range(list_spheres.shape[0]):
    axs.add_patch(plt.Circle((list_spheres[i, 0], list_spheres[i, 1]), 0.5e-3, color='r', fill=False))
axs.set_xlim([-0.002, 0.006])
axs.set_ylim([-0.003, 0.006])
plt.show()