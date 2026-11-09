import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.auto import tqdm

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[1]))  # repository root
from utils import *  # key functions for this project
from paths import ANALYSIS_DATA_FOLDER  # data locations; see paths.py

THRESHOLD = 50

print('Getting distributions for looped/unlooped separately')
loop_nums_CTCF = np.array([0, 1, 2])
noise_levels = np.array([0, 10, 20, 30, 40, 50])
reps = np.arange(10,19+1)

rows = []  # counts of occurrences of each nonnegative integer: 0, 1, 2 nm...

total_iters = len(loop_nums_CTCF) * len(noise_levels) * len(reps)

with tqdm(total=total_iters, desc="gathering 3D distances", unit="rep") as pbar:
    for loop_num in loop_nums_CTCF:
        for noise in noise_levels:
            
            dist_3D_looped = []
            dist_3D_unlooped = []
                
            for rep in reps:              
                traj = load_trajectory(loop_num, rep, noise)
                ctcf_state = load_ctcf_state(loop_num, rep)
                dist_3D_looped += list(traj[ctcf_state==1])
                dist_3D_unlooped += list(traj[ctcf_state==0])
                
                pbar.update(1)
                
            hist_looped = np.bincount(np.round(np.array(dist_3D_looped)).astype('int'))
            hist_unlooped = np.bincount(np.round(np.array(dist_3D_unlooped)).astype('int'))
                
            rows.append(
                {
                    "loop": loop_num,
                    "noise": noise,
                    "hist_looped": hist_looped,
                    "hist_unlooped": hist_unlooped
                }
            )
            
dist_3D_distribution_histograms_with_looping = pd.DataFrame(rows)

dist_3D_distribution_histograms_with_looping.to_pickle(ANALYSIS_DATA_FOLDER / '3D_dist_histograms_with_looping.pkl')
