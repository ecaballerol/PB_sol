import os
os.environ["OMP_NUM_THREADS"] = "7" # export OMP_NUM_THREADS=4
import numpy as np
from os.path import join,exists,basename
import math
from glob import glob

#Define the UTM
utmzone = 19

#Set Preliminar Hypocenter
lat_hypo = -31.637
lon_hypo = -71.741
dep_hypo =  23.3
M0_est = 8.3
#+++++++++++++++++++++++++++++++++++
#DATA dir
gps_dir = './DATA/GNSS'
enu = True
exclude_distance = None
gps_factor = 1 #1e-3 in case of gps in mmm

insar_dir = './DATA/INSAR'


#+++++++++++++++++++++++++
# Fault Geometry Parameters
GeometryFile = False
fault_dir = '../GEOMETRY'

strike = 4
dip = 22
FaultGeo = {'strike':strike,'dip':dip,\
            'n_strike':10,'n_dip':5,\
            'grid_size':2}

TopEdge = {'depth':5} 

comp_GFs    = True
GFdir       = './GFs' # Green's function directory
VelModel = 'okada'

plot_figs   = True
