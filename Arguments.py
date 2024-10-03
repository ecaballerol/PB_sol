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
#+++++++++++++++++++++++++++++++++++++
#DATA dir
sta_data = False
gps_dir = './DATA/GNSS'
enu = True
exclude_distance = None
gps_factor = 1 #1e-3 in case of gps in mmm

insar_dir = './DATA/INSAR'

# HR GPS data
seis_dir_hrgps       = './DATA/HRGPS'
sac_lst_hrgps        = join(seis_dir_hrgps,'o_inv_sac_file_lst')

# Strong motion data
seis_dir_sm       = './DATA/SEISMIC'
sac_lst_sm        = join(seis_dir_sm,'o_inv_sac_file_lst')


#+++++++++++++++++++++++++++++++++++++++
# Fault Geometry Parameters
GeometryFile = False
fault_dir = '../GEOMETRY'

strike = 4
dip = 22
FaultGeo = {'strike':strike,'dip':dip,\
            'n_strike':4,'n_dip':4,\
            'grid_size':4}


comp_GFs    = False #If true, will erase everythin in dir
comp_KinGFs = False
comp_bigG = True
GFdir       = './GFs' # Green's function directory
VelModel = 'okada'

plot_figs   = True

QuickInv = True
TypeInv = 'constrain'

# Green's functions for seismic data
earth_model = '/Users/EmmanuelC/Documents/Wphase_DB/valparaiso.model'
gf_db_dir   = '/Users/EmmanuelC/Documents/Wphase_DB/GF_JAPAN00'
gf_db_scale = 1.0e-28 # 1.0e-26 in cm -> 1.0e-28 in m
Ntriangles  = 40
Dtriangles  = 1.
Npt         = 10
Nmesh       = 10
Vr = 1.75 #km/s
# Kernel functions
seis_ref_model  = 'chile'
kernel_db_dir   = '../kernels/kinematic/Kernel_Pisagua_2014'
kernel_db_scale = 1.0e-28 # 1.0e-26 in cm -> 1.0e-28 in m

# Filter parameters for seismic data
#BP    = [0.009,0.05]
BP    = [0.01,0.06667] #v1 data
ORDER = 4
