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

def calcDim(Mest):
    '''
    Function to calculate the Area and width of an earthquake based on Wells et al. (1994)
    length = extension along strike
    width = extension along dip

    '''
    logStr = -2.44 + Mest * 0.59
    logRA = -3.49 + Mest * 0.91
    logDip = -1.01 + M0_est * 0.32
    #We add 20 percent 
    length = 1.2 * (10 ** logStr )
    area = 1.2 * (10 ** logRA )
    width = 1.6 * (10 ** logDip )
    return area, length, width

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
    if FaultGeo['strike'] >= 360:
        FaultGeo['strike'] = FaultGeo['strike'] -  360
    area,length,width = calcDim(M0_est)
    FaultGeo['length']=length
    FaultGeo['width'] = width
    #if FaultGeo['strike'] >0 and FaultGeo['strike'] < 90:
    Toplon = lon_hypo - (FaultGeo['width']/2 *np.cos(np.deg2rad(strike)))/111 
    Toplat = lat_hypo + (FaultGeo['width']/2 *np.sin(np.deg2rad(strike)))/111 
    #elif FaultGeo['strike']>=90 and FaultGeo['strike']<180:
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
            if enu:
                GPS.read_from_enu(GPSfile,header=1,factor=gps_factor)
            else:
                GPS.read_from_en(GPSfile,header=1,factor=gps_factor)

            if exclude_distance is not None:
                print('Excluding stations further than: ' + str(exclude_distance))
                GPS.reject_stations_awayfault(exclude_distance, fault)
            if enu:
                GPS.buildCd(direction='enu')
            else:
                GPS.buildCd(direction='en')
            data_avail.extend([GPS])
            count +=1
else:
    print('not GNSS available')

#InSAR directory
if os.path.isdir(insar_dir):
    if os.listdir(insar_dir):
        count =0
        for ifile in os.listdir(insar_dir):
            if ifile.endswith('.rsp'):
                print('INSAR available')
                InSAR = ir('INSAR, dataset: ' + str(count),utmzone=utmzone)
                InSARfile= os.path.join(insar_dir,os.path.splitext(ifile)[0])
                InSAR.read_from_varres(InSARfile,factor=0.01,step=0,header=2,cov=False)
                InSAR.buildCd(sigma=0.00605,lam=7.750)
                data_avail.extend([InSAR])
                count +=1
    else:
        print('not file found')
else:
    print('not InSAR available')



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
        elif data.dtype == 'gps':            
            fault.buildGFs(data,method=VelModel,vertical=enu)
        else:
            fault.buildGFs(data,method=VelModel)

    else: # Load GFs
        GfSS = os.path.join(GFdir,'{}_{}_{}.gf'.format(fault.name, data.name, 'SS'))
        GfDS = os.path.join(GFdir,'{}_{}_{}.gf'.format(fault.name, data.name, 'DS'))
        if data.dtype == 'gps':
            fault.setGFsFromFile(data,strikeslip=GfSS,dipslip=GfDS,vertical=enu)
        elif data.dtype == 'tsunami':
            G_SS,G_DS = dart.getGF(dart_wav,fault,factor=0.01)
            fault.setGFs(data,strikeslip=[G_SS],dipslip=[G_DS])
        else:
            fault.setGFsFromFile(data,strikeslip=GfSS,dipslip=GfDS,vertical=True)

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
    if dataset.dtype == 'gps':
        dataset.buildsynth(slv.faults,vertical=enu)
    else:
        dataset.buildsynth(slv.faults)

# %%
# Plot arguments
if plot_figs:
    plt.close('all')
    for idata in data_avail:
    # Prepare figures
        gp = geoplt(lonmin=lon_hypo-1.5,lonmax=lon_hypo+1.5,latmin=lat_hypo-2.,latmax=lat_hypo+2.,figsize=[(8,9),(8,9)])
    # Plot dipslip
        gp.faultpatches(fault, slip='dipslip', colorbar=True,plot_on_2d=True,alpha=0.6,cmap='Reds',cbaxis=[0.7, 0.4, 0.1, 0.02])
        gp.setzaxis(depth=int(fault.depth*1.3),zticklabels=None)
    # Plot gps data/predictions
        if idata.dtype == 'gps':
            gp.gps(idata,scale=5.,legendscale=1,data=['data','synth'],color=['k','b'])
            gp.gpsverticals(idata,data=['data','synth'],markersize=[100,30],norm=[-0.3,0.3],cbaxis=[0.7, 0.2, 0.1, 0.02])
        elif idata.dtype=='insar':
            gp.insar(idata,data='res',plotType='decimate',colorbar=True,alpha=0.7)
        gp.carte.coastlines()
        gp.carte.plot(fault.hypo_lon,fault.hypo_lat,'*k',ms=13,zorder=5)
        title= idata.name + ' fit'
        gp.titlemap(title)
        gp.faille.set_title('Slip')
    gp.show()

# %%
