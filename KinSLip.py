#%%

# Rapid Estimation of Slip 

#Created by E.Caballero GNS

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import shutil
import scipy.signal as signal
import copy

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

## Import Wavemod
import wavemod as wm

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
    logDip = -1.01 + Mest * 0.32
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
                    FaultGeo['grid_size'],FaultGeo['n_strike'],FaultGeo['n_dip'],leading='dip')
    Nstrike = fault.f_nstrike
    Ndip = fault.f_ndip
fault.setTrace(delta_depth=0.5)
fault.computeArea()
fault.trace2ll()

# Set hypocenter
fault.setHypoXY(lon_hypo,lat_hypo,UTM=False) # Hypocenter (for preliminar solution.)

#See if we define the Mu for EDKS

print('NUMPATCH : ', len(fault.patch))

# Set mapping
fault.setFaultMap(Nstrike,Ndip,leading='dip',check_depth=True)

#%%
#Searching and loading Available Static DATA
data_avail = []
#We start with GNSS
if sta_data==True:
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


#------------------------------------------
# SEISMIC DATA


seismic_data=[]
seis_std     = []
print("Init waveform engine")
gps_sac_files = []
gps_seis_std  = []
if os.listdir(seis_dir_hrgps):
    with open(sac_lst_hrgps) as f:
        for l in f:
            items = l.strip().split()        
            gps_sac_files.append(join(seis_dir_gps,items[0]))
            stdv = float(items[1])
            gps_seis_std.append(stdv)
    hrgps_data = seis('HR_GPS',utmzone=utmzone)
    hrgps_data.readSac(gps_sac_files)
    seismic_data.append(hrgps_data)
    seis_std.append(gps_seis_std)

sm_sac_files = []
sm_seis_std  = []
if os.listdir(seis_dir_sm):        
    with open(sac_lst_sm) as f:
        for l in f:
            items = l.strip().split()        
            sm_sac_files.append(join(seis_dir_sm,items[0]))
            stdv = float(items[1])
            sm_seis_std.append(stdv)
    sm_data = seis('Strong_motion',utmzone=utmzone)
    sm_data.readSac(sm_sac_files)
    
    seismic_data.append(sm_data)
    seis_std.append(sm_seis_std)

gauss_cor_std =[10.         ,6.]
velocity     = [False       ,False   ]

seismic_GFs=[]
for i in range(len(seismic_data)):
    data = seismic_data[i]    
    station_list = []
    station_lat = []
    station_lon = []
    for ifile in data.d:
        sacfile =data.d[ifile]
        if sacfile.kstnm not in station_list:
            print(sacfile.kstnm)
            station_list.append(sacfile.kstnm)
            station_lat.append(sacfile.stla)
            station_lon.append(sacfile.stlo)
    sm_data = seis('Strong_motion',utmzone=utmzone)
    sm_data.setStat(station_list,station_lon,station_lat)
    seismic_GFs.append(sm_data)

# %%
# Creating GFs directory
# Cleaning up

if comp_GFs == True: 
    if os.path.exists(GFdir):
        shutil.rmtree(GFdir)
    os.mkdir(GFdir)

#Create the GFs for each dataset
if sta_data==True:
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


#Create Kinematic GFs
if comp_KinGFs:
    for r in [0.,90.]:
        o_dir = os.path.join(GFdir,'rake_%.1f'%(r))
        if os.path.exists(o_dir):
            shutil.rmtree(o_dir)

faultbigG = copy.deepcopy(fault)
# Compute/Load Kinematic Green's functions                
for i in range(len(seismic_data)):
    data = seismic_data[i]    
    print("Compute Kinematic GFs for %s"%(data.name))    
    # Filter coefficients
    delta = data.d[data.sta_name[0]].delta
    sosfilter = signal.butter(ORDER,np.array(BP)*2*delta,'bandpass',output='sos')
    # Compute/Load/Save GFs
    if comp_KinGFs: 
        if GeometryFile == True:
        # Compute GFs
            # Waveform engine
            waveDB = wm.WaveDB(gf_db_dir,gf_db_scale)
            print(data.name)
            fault.buildKinGFsFromDB(data,waveDB,1.,0.,  filter_coef=sosfilter, differentiate=velocity[i])
            fault.buildKinGFsFromDB(data,waveDB,1.,90., filter_coef=sosfilter, differentiate=velocity[i])
            # Save GFs
            fault.saveKinGF(data,outputDir=GFdir,rmdir=False)
        else:
            data=seismic_GFs[i]
            waveInt=wm.WaveInt(earth_model,npts=1024,delta=1.,T0=-20.)
            data.initWaveInt(waveInt)
            print(data.name)
            fault.buildKinGFs(data,35e9,0.,1.,rise_time=1,stf_type='dirac', out_type='V',filter_coef=sosfilter)
            fault.buildKinGFs(data,35e9,90.,1.,rise_time=1,stf_type='dirac', out_type='V',filter_coef=sosfilter)
            # Save GFs
            fault.saveKinGFs(data,o_dir=GFdir)
  
# %%
#Windowing the GFs database according to the seismic data
fault = copy.deepcopy(faultbigG)
for idat in range(len(seismic_data)):
    data=seismic_data[idat]
    try:
        if GeometryFile == True:
            fault.loadKinGF(data,[0.,90.],inputDir=GFdir)
        else:
            # data=seismic_GFs[i]
            fault.loadKinGFs(data,[0.,90.],i_dir=GFdir,components='individual')
    except:
        print('Not able to load kinematic GFs')
    fault.GFwindow(data)

#Creating Vr for each patch
fault.vr = np.zeros((len(fault.patch,)))
fault.vr[:]=3.5

# Build Big G matrix
print("Assemble bigG and bigD")
bigDfile = os.path.join(GFdir,'kinematicG.data.bin')
bigGfile = os.path.join(GFdir,'kinematicG.gf.bin')
fault.setBigDmap(data)

if comp_bigG:
    #Add FastSweep
    fastSweep = wm.FastSweep()
    tmin =fault.buildBigGD(fastSweep,seismic_data,[0.,90.],3.7,Ntriangles,Dtriangles,dtype='float32',\
                           fastsweep=True,indexing='Altar')
#Create the bigG for the kinematic problem

#Create the Cd
for i in range(len(seismic_data)):
    data = seismic_data[i]
    data.buildDiagCd(seis_std[i])
fault.buildBigCd(seismic_data)
#+++++++++++++++++++++++++++++++++++++++++++++
#Quick inversion for big G
# Define the cost function
def costFunction(m, G, d, mprior,iCd=None, verbose=False):
    """
    Compute data + prior misfits.
    """
    dataMisfit = d - np.dot(G,m)
    if iCd is None:
        dataLikely = np.dot(dataMisfit, dataMisfit.T)
    else:
        dataLikely = np.dot(dataMisfit, np.dot(iCd,dataMisfit.T))
    priorMisfit = m - mprior
    priorLikely = np.dot(priorMisfit, priorMisfit.T)
    if verbose:
        print(0.5 * dataLikely)
#return 0.5 * dataLikely + 0.5 * priorLikely
    return 0.5 * dataLikely

if QuickInv:
    import scipy.linalg as scilin
    import scipy.optimize as sciopt
    from scipy.optimize import minimize, nnls
    import matplotlib.patches as mpl_patches
    
    D = fault.bigD
    G = fault.bigG
    Cd = data.Cd
    iCd = scilin.inv(Cd)

    Np = len(fault.patch)
    print("Quick LSQ inversion")
    if TypeInv == 'positivity':
        mpost, rnorm = sciopt.nnls(G, D)
        figopt = 'positivity_'
    elif TypeInv == 'norestrain':
        mpost = np.dot( np.dot( scilin.inv(np.dot( G.T, G )), G.T ), D)
        figopt = 'no_restrain_'
    elif TypeInv== 'constrain':
        bounds = []
        for it in range(Ntriangles):
            for ip in range (Np):
                bounds.append((-2,2))
            for ip in range (Np):
                bounds.append((0.0,35.0))
        maxfun = 10000
        iterations=500
        checkIter= True
        options = {'disp': checkIter, 'maxiter': iterations}
        options['maxfun']= maxfun
        method='L-BFGS-B'
        mprior = np.zeros((len(fault.patch)*Ntriangles*2,))
        checkNorm = True
        constraints = ()
        res = minimize(costFunction, mprior,args=(G,D,mprior,iCd,checkNorm),constraints=constraints,method=method, bounds=bounds,
               options=options)
        mpost = res.x
        figopt = 'constrain_'


#Synthetics
P= np.dot(G,mpost)

from scipy import integrate
tmin = np.around(np.array(tmin)).astype(int)

#Arranging m vector to Slip rate arrays
MrateStr = np.zeros((Np,Ntriangles))
MrateDip = np.zeros((Np,Ntriangles))
Ndip = fault.f_ndip
Nstrike = fault.f_nstrike

for it in np.arange(Ntriangles):
    for ip in np.arange(Np):
        istr = 2 * it * Np + ip
        idip = ((2 * it +1 )* Np ) + ip
        MrateStr[ip,it] = mpost[istr]
        MrateDip[ip,it] = mpost[idip]

SS = integrate.cumulative_trapezoid(MrateStr,initial=0,dx=Dtriangles,axis=1)
SD = integrate.cumulative_trapezoid(MrateDip,initial=0,dx=Dtriangles,axis=1)
SlipTot = np.sqrt(SS ** 2 + SD ** 2)

fault.slip[:,0]=SS[:,-1]
fault.slip[:,1]=SD[:,-1]

tmax = max(tmin) + Ntriangles +10
SStrike_Plane = np.zeros((Ndip,Nstrike,tmax))
SDip_Plane = np.zeros((Ndip,Nstrike,tmax))
Slip_Plane = np.zeros((Ndip,Nstrike,tmax))

STF = np.zeros(tmax) 
for ip in np.arange(Np):
    STF[tmin[ip]:tmin[ip]+Ntriangles] += np.sqrt(MrateDip [ip,:] **2 + MrateStr[ip,:] ** 2)

    #plotting the STF
    
mu=35e9    
STF = STF * mu * fault.area[0] * 1e6 

M0 = integrate.trapezoid(STF)
Mw = 2./3.*(np.log10(M0)-9.1)

handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white",lw=0, alpha=0)] * 2
plt.plot(STF)
plt.xlabel('time, s')
plt.ylabel('N m/s')
labels = []
labels.append("M0 = {0:.2e}".format(M0))
labels.append("Mw = {0:.3g}".format(Mw))
plt.legend(handles,labels, loc='best',fancybox=True, framealpha=0.7,handlelength=0, handletextpad=0)
plt.grid()
figurefile = './Figures/STF_' + figopt +  '.pdf'
plt.savefig(figurefile,bbox_inches='tight')
plt.close()

gp = geoplt(lonmin=lon_hypo-1.5,lonmax=lon_hypo+1.5,latmin=lat_hypo-2.,latmax=lat_hypo+2.,figsize=[(8,9),(8,9)]) 
gp.faultpatches(fault, slip='total', colorbar=True,plot_on_2d=True,alpha=0.6,cmap='Reds',cbaxis=[0.7, 0.4, 0.1,0.02])
gp.setzaxis(depth=int(fault.depth*1.3),zticklabels=None)
gp.carte.plot(fault.hypo_lon,fault.hypo_lat,'*k',ms=13,zorder=5)
gp.carte.coastlines()


'''
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
'''