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
    count =0
    for ifile in os.listdir(gps_dir):
        if ifile.endswith('.dat'):
            print('GNSS available')
            GPS = gr('GPS, dataset: ' + str(count),utmzone=utmzone)
            GPSfile= os.path.join(gps_dir,ifile)
            GPS.read_from_enu(GPSfile,header=1)
            GPS.buildCd(direction='enu')
            data_avail.extend([GPS])
            count +=1
else:
    print('not GNSS available')

# %%
#Initialize fault object
# Two approaches, either Initialize planar fault, either 
#look for the geometry
if GeometryFile is True:
    for ifile in os.listdir(fault_dir):
        # Initialize a fault from txt
        fault = rectangular_fault('EarSlinv', utmzone=utmzone)

        # Get the fault patches from gocad
        fault.read3DsquareGrid(os.path.join(fault_dir,ifile))
        
        #fault.initializeslip()
        fault.initializekinmodel()
else:
    print('Creating a planar fault')
    if FaultGeo['strike'] >0 and FaultGeo['strike'] < 90:
        Toplon = lon_hypo - (FaultGeo['width']/2 *np.cos(np.deg2rad(strike)))/111 
        Toplat = lat_hypo + (FaultGeo['width']/2 *np.sin(np.deg2rad(strike)))/111 
    #elif FaultGeo['strike']>=90 and FaultGeo['strike']<180:
#    Toplon = lon_hypo - (FaultGeo['width']/2 *np.cos(np.deg2rad(strike)))/111 
#    Toplat = lat_hypo + (FaultGeo['width']/2 *np.sin(np.deg2rad(strike)))/111 
    TopEdge['lon']=Toplon
    TopEdge['lat']=Toplat
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

# %%
#Assemble the inverse problem

#  Assemble GFs
# fault.assembleGFs(static_datasets,slipdir='sd',polys=poly)
fault.assembleGFs(data_avail,slipdir='sd',polys=None)

# Assemble the datasets
fault.assembled(data_avail)
fault.assembleCd(data_avail)

fault.buildCm(sigma=10.0,lam= 20.0) #See Radiguet et al., 2010
#fault.buildCm(sigma=10.0,lam= 30.0) #See Radiguet et al., 2010

# Assemble the problem
slv = multiflt('FastInv', [fault])
slv.assembleGFs()
slv.assembleCm()

# Write matrices to binary files for Altar
slv.writeGFs2BinaryFile(outfile=os.path.join(GFdir,'static.gf'),dtype='d')
slv.writeData2BinaryFile(outfile=os.path.join(GFdir,'static.data'),dtype='d')
slv.writeCd2BinaryFile(outfile=os.path.join(GFdir,'static.Cd'),dtype='d', scale=1.0)
slv.writePatchAreasFile(outfile=os.path.join(GFdir,'variableDip_patchAreas.dat'))

# %%
#First Slip Inversion

# Slip boundaries
bounds = []
for i in range(fault.N_slip): # Along strike
    bounds.append((-2.,2.))
for i in range(fault.N_slip): # Along dip
    bounds.append((0.,30.))

# Inversion
slv.ConstrainedLeastSquareSoln(bounds=bounds,method='L-BFGS-B',checkIter=True)
# Distribute slip values in fault
slv.distributem()
# Write model in a text file
np.savetxt('static_LSQ.txt',slv.mpost)    


# Compute predictions
for dataset in data_avail:
    dataset.buildsynth(slv.faults)

# %%
# Plot arguments
if plot_figs:
    plt.close('all')
    for idata in data_avail:
    # Prepare figures
        gp = geoplt(lonmin=-73,lonmax=-69,latmin=-33,latmax=-29,figsize=[(8,9),(8,9)])
    # Plot dipslip
        gp.faultpatches(fault, slip='dipslip', colorbar=True,plot_on_2d=True,alpha=0.6,cmap='Reds',cbaxis=[0.7, 0.4, 0.1, 0.02])
        gp.setzaxis(depth=50,zticklabels=None)
    # Plot gps data/predictions
        if idata.dtype == 'gps':
            gp.gps(idata,scale=5.,legendscale=1,data=['data','synth'],color=['k','b'])
            gp.gpsverticals(idata,data=['data','synth'],markersize=[100,30],norm=[-0.3,0.3],cbaxis=[0.7, 0.2, 0.1, 0.02])
        gp.carte.coastlines()
        gp.carte.plot(fault.hypo_lon,fault.hypo_lat,'*k',ms=13,zorder=5)
        title= idata.name + ' fit'
        gp.titlemap(title)
        gp.faille.set_title('Slip')
    gp.show()

# %%
