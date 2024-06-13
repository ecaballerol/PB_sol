import os
os.environ["OMP_NUM_THREADS"] = "7" # export OMP_NUM_THREADS=4
import numpy as np
from os.path import join,exists,basename
import math
from glob import glob


#DATA dir

gps_dir = '../DATA/GNSS'