import os
from os.path import join,exists,basename

from glob import glob


# Import externals
import matplotlib
matplotlib.use('PDF')
import sys
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import scipy.signal as signal
from scipy import interpolate
from scipy import integrate
import h5py
# Import personal libraries
import csi.RectangularPatchesKin as rectangular_fault # Use only for kinematic modeling
import csi.multifaultsolve as multiflt
import csi.faultpostproc as faultpostproc
import csi.geodeticplot as geoplt
import csi.insar as ir
import csi.gps   as gr
import csi.tsunami    as tsunami
import csi.seismic    as seis

## Import Wavemod
import wavemod as wm

#!/usr/bin/env python
LOGDIR        = 'LOG'
DLAT, DLON    = 4., 4.        # Half-size of the map region
OPDFFILE      = 'inversion_pages.pdf'

FIGSIZE   = [11.69,8.270]
#FIGSIZE   = [5.84,4.135]
YLIM_AUTO = True
YLIMFIXED = [-9,12] # Y lim if YLIM_AUTO = False
NC = 3 # Number of columns
NL = 5 # Number of lines


def show_basemap(ax,evla,evlo,stla,stlo,coords,flagreg=True):
    
    if flagreg:
        projection = ccrs.PlateCarree()
    else:
        projection = ccrs.cartopy.crs.Orthographic(central_longitude=evlo, central_latitude=evla)
    #projection = ccrs.PlateCarree()
    
    pos  = ax.get_position().get_points()
    W  = pos[1][0]-pos[0][0] ; H  = pos[1][1]-pos[0][1] ;
    ax2 = plt.axes([pos[1][0]-W*0.38,pos[0][1]+H*0.01,H*1.08,H*1.00],projection=projection)
    
    if flagreg:
        ax2.set_extent([evlo-DLON,evlo+DLON,evla-DLAT,evla+DLAT],crs=projection)
    else:
        radius_m = 15 * 111e3  # 1 degree ~ 111 km
        
        # Set the limits for the Orthographic projection
        ax2.set_xlim(-radius_m, radius_m)
        ax2.set_ylim(-radius_m, radius_m)
        #ax2.set_global()
    
    #ax2.set_extent([evlo-DLON,evlo+DLON, evla-DLAT,evla+DLAT],crs=projection)
    
    ax2.coastlines(linewidth=0.5)
    ax2.add_feature(cfeature.LAND,facecolor='0.75')
    if flagreg:
        gl = ax2.gridlines(draw_labels=True,crs=ccrs.PlateCarree())
    else:
        gl = ax2.gridlines(crs=ccrs.PlateCarree())
    gl.top_labels = False
    gl.left_labels = False
    gl.bottom_labels = False
    gl.right_labels = False
    ax2.plot(coords[:,1],coords[:,0],'o',color=(1.00000,  0.74706,  0.00000),ms=4.0,alpha=1.0,zorder=1000,transform=ccrs.PlateCarree())
    ax2.plot(stlo,stla,'o',color=(1,.27,0),ms=8,alpha=1.0,zorder=1001,transform=ccrs.PlateCarree())
    ax2.scatter(evlo,evla,c='b',marker=(5,1,0),s=120,zorder=1002,transform=ccrs.PlateCarree())
    
    return
    
    
plotparams2 = {'backend': 'pdf', 'axes.labelsize': 12, 'font.size': 12,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
        'legend.fontsize': 12, 'lines.markersize': 6, 'font.size': 12,  'savefig.dpi': 200,
            'keymap.back': ['left', 'c', 'backspace'], 'keymap.forward': ['right', 'v'],
            'keymap.fullscreen': 'f', 'keymap.grid': 'g', 'keymap.home': ['h', 'r', 'home'], 'keymap.pan': 'p',
            'keymap.save': 's', 'keymap.xscale': ['k', 'L'], 'keymap.yscale': 'l', 'keymap.zoom': 'o',
            'path.snap': True, 'savefig.format': 'pdf', 'pdf.compression': 9, 'figure.figsize': FIGSIZE}
plt.rcParams.update(plotparams2)

def data_plot(seismic_data,synthetic):
    cpt = ['red','gray']
    DLAT, DLON    = 4., 4.        # Half-size of the map region
    OPDFFILE      = 'inversion_pages.pdf'

    FIGSIZE   = [11.69,8.270]
    YLIM_AUTO = True
    YLIMFIXED = [-9,12] # Y lim if YLIM_AUTO = False
    NC = 3 # Number of columns
    NL = 5 # Number of lines
    ##PLOT IN WPHASE WAY
    nc = NC
    nl = NL
    title = 'kinematic data'
    perpage = nl*nc
    statnum = 0
    latstat = []
    lonstat = []
    for idata in seismic_data:
        statnum += len(idata.sta_name)
        latstat.extend(idata.lat)
        lonstat.extend(idata.lon)
    latstat =np.array(latstat)
    lonstat = np.array(lonstat)
    ntot   = statnum
    HRstat = len(seismic_data[0].sta_name)
    coords = []
    coords = np.array([latstat,lonstat]).T
    npages = np.ceil(float(ntot)/float(perpage))
    nchan = 1
    count = 1
    pages = 1
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.08,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.55)
    print ('%d pages:'%(npages))
    pp = matplotlib.backends.backend_pdf.PdfPages(OPDFFILE)
    m = None


    di = 0
    for l in np.arange(statnum):
        if l < HRstat:
            data = seismic_data[0]
            idx = l
        else:
            data = seismic_data[1]
            idx = l - HRstat
        chan = data.d[data.sta_name[idx]].kcmpnm
        loc  = data.d[data.sta_name[idx]].khole


        if count > perpage:
    #        plt.suptitle(title+ ',   p %d/%d'%(pages,npages), fontsize=16, y=0.95)
            ofic = 'page_W_%02d.pdf'%(pages)
            print (ofic)
            fig.set_rasterized(True)
            pp.savefig(orientation='landscape')#,format='pdf')
            plt.close()
            pages += 1
            count = 1
            fig = plt.figure()
            fig.subplots_adjust(bottom=0.08,top=0.87,left=0.06,right=0.95,wspace=0.25,hspace=0.55)

        # Plot trace
        sac = data.d[data.sta_name[idx]]
        t1 = np.arange(sac.npts,dtype='double')*sac.delta + sac.b - sac.o
        ax = plt.subplot(nl,nc,count)
        npts = sac.npts
        i = di
        l = npts + di
        plt.plot(t1,sac.depvar*1000.,'k',zorder=2)
        
        #Plot synthetics
        tmptr = synthetic[i:l]
        plt.plot(t1,tmptr*1000,'-',color=cpt[0],alpha=0.8,zorder=1)

        plt.xlim([t1[0],t1[-1]+100])
        if YLIM_AUTO:
            a    = sac.depvar.min()*1000.
            b    = sac.depvar.max()*1000.
            ymin =  1.1*a
            ymax =  1.1*b
            if ymin>-1. or ymax<1.:
                ymin = -2.
                ymax =  2.
            ylims = [ymin,ymax]
            
        else:
            ylims = YLIMFIXED
        plt.ylim(ylims)
            
        # Annotations
        if sac.kcmpnm[2] == 'Z':
            label = r'%s %s %s %s $(\phi,\Delta) = %6.1f^{\circ}, %6.1f^{\circ}$'%(sac.knetwk,sac.kstnm, sac.kcmpnm, sac.khole,                                                                        sac.az, sac.gcarc)
        else:
            label  = r'%s %s %s %s $(\phi,\Delta,\alpha) = %6.1f^{\circ},'
            label += '%6.1f^{\circ}, %6.1f^{\circ}$'
            label  = label%(sac.knetwk,sac.kstnm, sac.kcmpnm, sac.khole,
                                sac.az, sac.gcarc, sac.cmpaz)
        plt.title(label,fontsize=12.0,va='center',ha='center')
        if not (count-1)%nc:
            plt.ylabel('mm',fontsize=12)
        if (count-1)/nc >= nl-1 or nchan+nc > ntot:
            plt.xlabel('time, sec',fontsize=12)
        plt.grid()
        try:
            m = show_basemap(ax,sac.evla,sac.evlo,sac.stla,sac.stlo,coords,flagreg=False)
            pass
        except:
            #show_polarmap(ax,sac.az,sac.dist,coords)
            print ('No basemap module')
        count += 1
        nchan += 1
        di += npts
    ofic = 'page_W_%02d.pdf'%(pages)
    print (ofic)
    fig.set_rasterized(True)
    plt.suptitle(title + ',    p %d/%d'%(pages,npages), fontsize=16, y=0.95)
    pp.savefig(orientation='landscape')#,format='pdf')
    plt.close()
    pp.close()
    
    return
