"""
Locations of the data used by this repository.
Please download the data from Zenodo and place it in data/; it should include:

    simulation_data_raw/
    simulation_data_processed/
    analysis_data/
    2022_Fbn2_full_data.h5

By default the data folder is expected at `data/` in the root of this repository
(a folder or a symbolic link). To keep it elsewhere, either set the environment
variable LOOPING_THRESHOLDING_ANALYSIS_DATA to its location, or edit the default
path in DATA_ROOT below.
"""

import os
from pathlib import Path

# root folder of this repository
REPO_ROOT = Path(__file__).resolve().parent

# root folder of the data (see above)
DATA_ROOT = Path(os.environ.get("LOOPING_THRESHOLDING_ANALYSIS_DATA", REPO_ROOT / "data"))

SIM_RAW_FOLDER = DATA_ROOT / "simulation_data_raw"
SIM_PROCESSED_FOLDER = DATA_ROOT / "simulation_data_processed"
ANALYSIS_DATA_FOLDER = DATA_ROOT / "analysis_data"
FBN2_DATA_FILE = DATA_ROOT / "2022_Fbn2_full_data.h5"

# table of the simulated loops (anchor monomers, loop type), part of this repository
LOOP_INFO_FILE = REPO_ROOT / "loop_info.csv"
