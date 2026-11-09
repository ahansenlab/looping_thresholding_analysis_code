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

loop_nums = np.array([0, 1, 2, 3, 4, 5])
looping_probs = pd.Series(index=loop_nums, dtype='float', name='looping_prob')
looping_probs.index.name = 'loop_num'
n_reps = 19 - 10 + 1  # 10..19 inclusive

## -- CTCF loops: as defined by cohesin positions --

loop_nums_CTCF = np.array([0, 1, 2])
total_iters = len(loop_nums_CTCF) * n_reps

print('Working on CTCF loops')
with tqdm(total=total_iters, desc="looping probabilities (CTCF)", unit="rep") as pbar:
    for loop_num in loop_nums_CTCF:
        looping_probs.loc[loop_num] = prob_loop_CTCF_all_reps(loop_num, pbar=pbar)
        
## -- EP loops: prob. under 50nm minus background --

loop_nums_EP = np.array([3, 4, 5])
total_iters = len(loop_nums_EP) * n_reps * 2  # extra factor of 2 for sticky vs. non-sticky

print('Working on EP loops')
with tqdm(total=total_iters, desc="looping probabilities (EP)", unit="rep") as pbar:
    for loop_num in loop_nums_EP:
        
        prob_under_threshold_regular = prob_under_threshold_all_reps(loop_num, noise=0, non_sticky=False, threshold=THRESHOLD, pbar=pbar)
        prob_under_threshold_non_sticky = prob_under_threshold_all_reps(loop_num, noise=0, non_sticky=True, threshold=THRESHOLD, pbar=pbar)
        
        looping_prob_this_loop = prob_under_threshold_regular - prob_under_threshold_non_sticky
        
        looping_probs.loc[loop_num] = looping_prob_this_loop

looping_probs.to_csv(ANALYSIS_DATA_FOLDER / 'looping_probs.csv', float_format="%.6f")