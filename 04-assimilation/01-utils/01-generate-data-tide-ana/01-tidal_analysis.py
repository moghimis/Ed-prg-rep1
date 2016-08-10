#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
script to read water elevation data and make them ready for tidal analysis

"""
from ctd.plotting import plot

__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2016, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import netCDF4
import netcdftime
from datetime import datetime,timedelta
from dateutil.parser import *
import scipy.io as sio
from   collections  import defaultdict
import numpy as np
import os,sys
import pandas as pd
import cPickle as pickle
import string
import matplotlib.pyplot as plt
import  pynmd.plotting.plot_settings as ps
from  pynmd.plotting.plot_settings import mat2py_datenum
import glob
from pynmd.plotting.plot_settings import jetWoGn


from pynmd.tools.tide_analysis import read_roms_zeta
from pynmd.tools.tide_analysis import do_r_t_tide_analysis
from pynmd.tools.tide_analysis import do_tappy_tide_analysis
from pynmd.tools.tide_analysis import projection_nri,lonlat2xy


base_dir_pri = '/data01/01-projects/05-nasa-altimeter/03-model-outputs/01-NASA-try02/02-prior/00_run/'
base_dir_tru = '/data01/01-projects/05-nasa-altimeter/03-model-outputs/01-NASA-try02/01-true/00_run/'
#file_name   = '/out/zeta_x_y05_x10.nc'
file_name    = '/out/zeta_x_y05_x05.nc'

pickle_name  = '/out/rt_tide/r_t_tide_out_txt_files/r_t_tide_pickle/rt_tide_data.pickle'

base_dirs = [base_dir_pri , base_dir_tru] 
names     = ['prior'      , 'true']


picks = []
bats  = []    
for name,base_dir in zip(names,base_dirs):
    #pass collection to tidal analysis routine
    #constits = ['M2','N2','S2','Q1','K1','K2','P1','O1','M4']
    constits  = ['M2','K1']
    hisfile = base_dir + file_name
    if True:
        # Generate collection
        flow_data = defaultdict(dict)
        #
        read_roms_zeta(flow_data = flow_data,name=name,roms_hisfile=hisfile)
        #
        out_dir = base_dir + '/out/rt_tide/'
        rt_out_dir = do_r_t_tide_analysis  (flow_data = flow_data,constits=constits ,out_dir = out_dir)
        #tp_out_dir = do_tappy_tide_analysis(flow_data = flow_data,constits=constits ,out_dir = out_dir)

    pick_name  = base_dir + pickle_name
    picks.append(  pickle.load(open( pick_name , "r" ))  )

    #put together tide data 
    print '  >  ROMS ZETA  > ;'
    ncvar   = netCDF4.Dataset(hisfile,'r').variables
    dates   = netCDF4.num2date(times=ncvar['ocean_time'][:],units=ncvar['ocean_time'].units,calendar='standard')
    x       = ncvar['x_rho'][:]
    y       = ncvar['y_rho'][:]
    bats.append(ncvar['h']  [:])
    proj    = projection_nri()
    lon,lat = lonlat2xy(proj,x,y,inverse=True)
    elev    = ncvar['zeta'][:]
    nj,ni   = elev[0,:].shape
    mask    = elev[0,:].mask



tide_vars = ['amp','amp_err','pha','pha_err']
tide_matrix = []

for pick in picks:
    for iconst in range(len(constits)):
        const = constits[iconst]
        tide_data = np.zeros((4,nj,ni))
        tmp = 0
        for sta_name in pick.keys():
            i  = pick[sta_name]['i']
            j  = pick[sta_name]['j']
            df = pick[sta_name]['r_t_tide']
            #print i,j,df.loc[iconst][0]
            if df.loc[iconst][0]==const:
                tide_data[0,j,i] =  df.loc[iconst][1]   
                tide_data[1,j,i] =  df.loc[iconst][2]   
                tide_data[2,j,i] =  df.loc[iconst][3]   
                tide_data[3,j,i] =  df.loc[iconst][4]   
            
            else:
                print ' Consts are not arranged well!!!!'
        
        tmp = dict(amp = tide_data[0,:],amp_err= tide_data[1,:],pha= tide_data[2,:] * np.pi / 180.0 ,pha_err = tide_data[3,:])
        if  constits[iconst]=='M2':
            M2 = tmp
        else:
            K1 = tmp
        
    tide_matrix.append(dict(M2=M2,K1=K1))
            




      
#plots
#sys.exit()
##########
dpi = 600
figs_dir = string.join(base_dir_pri.split('/')[:-3],'/')+'/figs/'
os.system('mkdir -p ' +figs_dir )
#########
for const in constits:
    print 'plotting > ' , const
    #
    real1 = tide_matrix[1][const]['amp'] * np.cos(tide_matrix[1][const]['pha'])
    real0 = tide_matrix[0][const]['amp'] * np.cos(tide_matrix[0][const]['pha'])
    
    imag1 = tide_matrix[1][const]['amp'] * np.sin(tide_matrix[1][const]['pha'])
    imag0 = tide_matrix[0][const]['amp'] * np.sin(tide_matrix[0][const]['pha'])
    
    data  = np.sqrt( (real1 - real0)**2 +(imag1 - imag0)**2 )
    
    #data = np.sqrt( (tide_matrix[1][const][param] - tide_matrix[0][const][param])**2 )
    data = np.ma.masked_array(data,mask)
    data = np.ma.masked_array(data,np.isnan(data))
    vmax = data.max()
    
    # depth contours to plot
    levels = np.arange(0,vmax,vmax/50)   
    # tricontourf plot of water depth with vectors on top
    plt.figure(3,figsize=(10,8))
    plt.clf()
    plt.subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
    img = plt.contourf(lon, lat, data , levels = levels ,shading='faceted',cmap=plt.cm.jet,extend='both')
    #plt.axis(ax)
    plt.gca().patch.set_facecolor('0.5')
    cbar = plt.colorbar(img, shrink=0.8) #orientation = 'horizontal'
    cbar.set_label('rmse', rotation=90, fontsize=14)
    #Q = quiver(lonf[idv],latf[idv],uf[idv],vf[idv],scale=15,width=0.005)
    #maxstr='%3.1f m/s' % maxvel
    #qk = quiverkey(Q,0.2,0.9,maxvel,maxstr,labelpos='W')
    #plt.title('RMSE of the true and prior, %s, %s  \n ' % (const, param));
    plt.title('RMSE of the true and prior, %s \n ' % (const));
    plt.xlabel('Lon')
    plt.ylabel('Lat')
    plt.savefig(figs_dir+'/'+names[1]+'-'+names[0]+'_'+const+'.png',dpi=dpi)


ylim = plt.gca().get_ylim()
xlim = plt.gca().get_xlim()


#plot h_dif
def read_path(limits):
    data_files = ['/data01/01-projects/05-nasa-altimeter/01-progs/01-data/03-Satelite/01-tracks/data/Visu_TP_Tracks_HiRes_GE_V2.kml',
                 '/data01/01-projects/05-nasa-altimeter/01-progs/01-data/03-Satelite/01-tracks/data/Visu_TP_Interlaced_Tracks_HiRes_GE_V2.kml']
    
    lonmin,lonmax,latmin,latmax = limits
    #
    lat_reg = []
    lon_reg = []
    for data_file in data_files:
        print data_file
        finfo = open(data_file)
        il = 0
        for line in finfo:
            if ',0' in line: 
                if '<Point>' in line:
                    break
                words = string.split(line, ',')
                lon = float(words[0])
                lat = float(words[1])
                if (lat > latmin) & (lat < latmax) & (lon > lonmin) & (lon < lonmax):
                    #print lon , lat
                    lon_reg.append(lon)
                    lat_reg.append(lat)
            il += 1         
        finfo.close()

    return np.array(lon_reg),np.array(lat_reg)





data = (bats[1]- bats[0])
data = np.ma.masked_array(data,mask)
data = np.ma.masked_array(data,np.isnan(data))
vmax = data.max()

# depth contours to plot
levels = np.arange(-vmax,vmax,vmax/50)   
# tricontourf plot of water depth with vectors on top
plt.figure(3,figsize=(10,8))
plt.clf()
plt.subplot(111,aspect=(1.0/np.cos(np.mean(lat)*np.pi/180.0)))
img = plt.contourf(lon, lat, data , levels = levels ,shading='faceted',cmap=jetWoGn(),extend='both')
#plt.axis(ax)
plt.gca().patch.set_facecolor('0.5')
cbar = plt.colorbar(img, shrink=0.8) #orientation = 'horizontal'
cbar.set_label('rmse', rotation=90, fontsize=14)

# extent of new grid
lonlim1 = -79.0 ,-75.0
latlim1 =  31.0 , 35.1
limits  = np.array([lonlim1,latlim1]).flatten()
lon_sat,lat_sat = read_path(limits)
plt.scatter(lon_sat,lat_sat,s=20,c='r')
#Q = quiver(lonf[idv],latf[idv],uf[idv],vf[idv],scale=15,width=0.005)
#maxstr='%3.1f m/s' % maxvel
#qk = quiverkey(Q,0.2,0.9,maxvel,maxstr,labelpos='W')
plt.title('(true - prior) Bathymetry') #%s, %s  \n ' % (const, param));
plt.xlabel('Lon')
plt.ylabel('Lat')

plt.xlim(xlim)
plt.ylim(ylim)

plt.savefig(figs_dir+'/'+'bathy_diff.png',dpi=dpi)


