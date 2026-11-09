import numpy as np
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))  # repository root
from paths import SIM_PROCESSED_FOLDER  # data locations; see paths.py

sim_timestep = 0.02
monomer_size_nm = 21
step = 1

noise_nm = 0
for loop_num in [0, 1, 2, 3, 4, 5]:
    for rep in range(10,20):
        mat = np.loadtxt(SIM_PROCESSED_FOLDER / f'loop_{loop_num}_{rep}_step_{step}.txt')

        # rescale time and distance to seconds and nm
        mat[:,0] = mat[:,0] * sim_timestep
        mat[:,1:7] = mat[:,1:7] * monomer_size_nm

        # calculate the 3D distance
        dist_3D = np.sqrt(np.sum(np.square(mat[:,1:4]-mat[:,4:7]), 1))
        mat_with_dist_3D = np.concatenate((mat, dist_3D.reshape((-1,1))), axis=1)

        np.savetxt(SIM_PROCESSED_FOLDER / f'loop_{loop_num}_{rep}_step_{step}_noise_{noise_nm}nm.txt', mat_with_dist_3D, fmt='%.2f')


for noise_nm in [10, 20, 30, 40, 50]:
    for loop_num in [0, 1, 2, 3, 4, 5]:
        for rep in range(10,20):
            mat = np.loadtxt(SIM_PROCESSED_FOLDER / f'loop_{loop_num}_{rep}_step_{step}.txt')

            # rescale time and distance to seconds and nm
            mat[:,0] = mat[:,0] * sim_timestep
            mat[:,1:7] = mat[:,1:7] * monomer_size_nm

            # add noise (2*sigma for z-direction)
            mat[:,1] = mat[:,1] + np.random.normal(0, noise_nm, mat.shape[0])  # x1
            mat[:,2] = mat[:,2] + np.random.normal(0, noise_nm, mat.shape[0])  # y1
            mat[:,3] = mat[:,3] + np.random.normal(0, noise_nm*2, mat.shape[0])  # z1
            mat[:,4] = mat[:,4] + np.random.normal(0, noise_nm, mat.shape[0])  # x2
            mat[:,5] = mat[:,5] + np.random.normal(0, noise_nm, mat.shape[0])  # y2
            mat[:,6] = mat[:,6] + np.random.normal(0, noise_nm*2, mat.shape[0])  # z2

            # calculate the 3D distance
            dist_3D = np.sqrt(np.sum(np.square(mat[:,1:4]-mat[:,4:7]), 1))
            mat_with_dist_3D = np.concatenate((mat, dist_3D.reshape((-1,1))), axis=1)

            np.savetxt(SIM_PROCESSED_FOLDER / f'loop_{loop_num}_{rep}_step_{step}_noise_{noise_nm}nm.txt', mat_with_dist_3D, fmt='%.2f')
