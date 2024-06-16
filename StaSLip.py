#%%

# Rapid Estimation of Slip 

#Created by E.Caballero GNS

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import shutil

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
#We start with GNSS
if os.listdir(gps_dir):
    for ifile in os.listdir(gps_dir):
        if ifile.endswith('.dat'):
            print('GNSS available')
            GPS = gr('GPS',utmzone=utmzone)
            GPSfile= os.path.join(gps_dir,ifile)
            GPS.read_from_enu(GPSfile,header=1)
            GPS.buildCd(direction='enu')
            data_avail.extend([GPS])
else:
    print('not GNSS available')

# %%
#Initialize fault object
# Two approaches, either Initialize planar fault, either 
#look for the geometry

# Initialize a planar fault
fault = planar_fault('EarSlInv', utmzone=utmzone)

fault.buildFault(TopEdge['lon'],TopEdge['lat'],TopEdge['depth'],\
                FaultGeo['strike'],FaultGeo['dip'],\
                FaultGeo['length'],FaultGeo['width'],\
                FaultGeo['grid_size'],FaultGeo['n_strike'],FaultGeo['n_dip'])

fault.setTrace(delta_depth=0.5)
fault.computeArea()
fault.trace2ll()

# Set hypocenter
fault.setHypoXY(lon_hypo,lat_hypo,UTM=False) # Hypocenter (for preliminar solution.)

#See if we define the Mu for EDKS

print('NUMPATCH : ', len(fault.patch))


# %%
# Creating GFs directory
# Cleaning up

if comp_GFs == True: 
    if os.path.exists(GFdir):
        shutil.rmtree(GFdir)
    os.mkdir(GFdir)

#Create the GFs for each dataset

for data in data_avail:  # Loop over datasets
    if comp_GFs == True: # Compute GFs        

        if data.dtype == 'tsunami':
            G_SS,G_DS = dart.getGF(dart_wav,fault,factor=0.01)
            fault.setGFs(data,strikeslip=[G_SS],dipslip=[G_DS])
        else:            
            fault.buildGFs(data,method=VelModel,vertical=True)

    else: # Load GFs
        GfSS = os.path.join(GFdir,'{}_{}_{}.gf'.format(fault.name, data.name, 'SS'))
        GfDS = os.path.join(GFdir,'{}_{}_{}.gf'.format(fault.name, data.name, 'DS'))
        if data.dtype == 'gps':
            fault.setGFsFromFile(data,strikeslip=GfSS,dipslip=GfDS,vertical=True)
        elif data.dtype == 'tsunami':
            G_SS,G_DS = dart.getGF(dart_wav,fault,factor=0.01)
            fault.setGFs(data,strikeslip=[G_SS],dipslip=[G_DS])
        else:
            fault.setGFsFromFile(data,strikeslip=GfSS,dipslip=GfDS,vertical=False)

# # Write Green's functions
if comp_GFs == True: # Write GFs
    fault.saveGFs(outputDir=GFdir)

#  Assemble GFs
# fault.assembleGFs(static_datasets,slipdir='sd',polys=poly)
fault.assembleGFs(data_avail,slipdir='sd',polys=None)




# %%
