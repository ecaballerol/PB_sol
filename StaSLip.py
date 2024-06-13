#%%

# Rapid Estimation of Slip 

#Created by E.Caballero GNS

import numpy as np
import matplotlib.pyplot as plt
import os


#Import CSI libraries
import csi.insar as ir
import csi.gps   as gr
import csi.tsunami    as tsunami
import csi.seismic    as seis

#------------------------------------------------------------------
#------------------------------------------------------------------
# Load Arguments
from Arguments import *

#%%
#Searching Available DATA

if os.listdir(gps_dir):
    print('GNS')
else:
    print('notGNS')

#Loading data

# %%
#Initialize fault object

