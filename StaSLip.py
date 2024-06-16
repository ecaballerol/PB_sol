#%%

# Rapid Estimation of Slip 

#Created by E.Caballero GNS

import numpy as np
import matplotlib.pyplot as plt
import os


#Import CSI libraries Data
import csi.insar as ir
import csi.gps   as gr
import csi.tsunami    as tsunami
import csi.seismic    as seis

#Import CSI libraries Fault
import csi.RectangularPatchesKin as rectangular_fault # Use only for kinematic modeling
import csi.multifaultsolve as multiflt
import csi.geodeticplot as geoplt
import csi.planarfaultkinematic as planar_fault

#------------------------------------------------------------------
#------------------------------------------------------------------
# Load Arguments
from Arguments import *

#%%
#Searching and loading Available DATA
data_avail = []
if os.listdir(gps_dir):
    print('GNS')
    GPS = gr('GPS',utmzone=utmzone)
    GPS.read_from_enu(gps_dir,header=1)#,factor=1000.)
    GPS.buildCd(direction='enu')
else:
    print('notGNS')

#Other approach
# for idir in os.listdir('DATA'):
#     print(idir)
#Loading data

# %%
#Initialize fault object
# Two approaches, either Initialize planar fault, either 
#look for the geometry
# Initialize a planar fault
fault = planar_fault('EarSlInv', utmzone=utmzone)


