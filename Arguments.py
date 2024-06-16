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

#DATA dir

gps_dir = '../DATA/GNSS'

print('tete')

#+++++++++++++++++++++++++
# Fault Geometry Parameters

TopEdge = {'lon': -72.419,'lat':-32.25,'depth':10} 
FaultGeo = {'length':300,'width':180,\
            'strike':4,'dip':22,\
            'n_strike':10,'n_dip':10,\
            'grid_size':2}


comp_GFs    = True
GFdir       = './GFs' # Green's function directory