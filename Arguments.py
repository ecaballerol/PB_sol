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

#+++++++++++++++++++++++++
# Fault Geometry Parameters
strike = 4
dip = 22
FaultGeo = {'length':320,'width':180,\
            'strike':strike,'dip':dip,\
            'n_strike':10,'n_dip':10,\
            'grid_size':2}
#if FaultGeo['strike'] >0 and FaultGeo['strike'] < 90:
Toplon = lon_hypo - (FaultGeo['width']/2 *np.cos(np.deg2rad(strike)))/111 
Toplat = lat_hypo + (FaultGeo['width']/2 *np.sin(np.deg2rad(strike)))/111 
#elif FaultGeo['strike']>=90 and FaultGeo['strike']<180:
#    Toplon = lon_hypo - (FaultGeo['width']/2 *np.cos(np.deg2rad(strike)))/111 
#    Toplat = lat_hypo + (FaultGeo['width']/2 *np.sin(np.deg2rad(strike)))/111 

TopEdge = {'lon': Toplon,'lat':Toplat,'depth':10} 

comp_GFs    = True
GFdir       = './GFs' # Green's function directory
VelModel = 'okada'

plot_figs   = True
