#%%

# Rapid Estimation of Slip 

#Created by E.Caballero GNS

import numpy as np
import matplotlib.pyplot as plt
import os

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

