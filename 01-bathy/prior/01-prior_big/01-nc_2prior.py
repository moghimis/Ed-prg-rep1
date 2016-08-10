#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
"""
script to read Gebco data and generate ROMS grid

"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2016, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import netCDF4
from numpy import *
import os,sys
import glob
import pylab as pl
import scipy.io as sio
from octant.grid import *
import octant.csa as csa

def interpg(x_old,y_old,data_old,x_new,y_new):
    csa_interp = csa.CSA(x_old,y_old,data_old)
    data_new = csa_interp(x_new,y_new)
    return data_new


def interp3(x,y,b,xnew,ynew,method):
    xf=x.flatten()
    yf=y.flatten()
    bf=b.flatten()
    interpm=method
    if interpm=='tri':
        xnewf=xnew.flatten()
        ynewf=ynew.flatten()
        print 'tri interp ...'
        from delaunay import  triangulate
        tri=triangulate.Triangulation(xf, yf)
        interp_b=tri.nn_extrapolator(bf)
        bnewf = interp_b(xnewf,ynewf)
        bnew=bnewf.reshape(xnew.shape)
    elif interpm=='csa' :
        print 'csa interp ...'
        import octant.csa as csa
        csa_interp = csa.CSA(xf, yf,bf)
        bnew = csa_interp(xnew,ynew)
    elif interpm=='grd' :
        print 'griddata interp ...'
        bnew=pl.griddata(xf,yf,bf,xnew,ynew)
    return bnew


#Cut to smaller region
i1=20
i2=130
#
j1=40
j2=140

print i1,i2,j1,j2

k=1

land=-10.0
#inp_prior  = 'inp/nri_coastal_grd_roms_v2etopo2_.nc'
inp_prior  = 'inp/nri_coastal_grd_roms_v2gebco_.nc'  
ncf_parent = 'inp/nri_coastal_grd_roms_v2gebco_.nc'  
new_prior  = 'new_prior.nc'

comm='cp  '+ncf_parent+'  tmp.nc'
os.system(comm)
 
################## topo 200 read in #############
tnc=netCDF4.Dataset(ncf_parent)
ncv=tnc.variables
xc  = ncv['x_rho'][:]
yc  = ncv['y_rho'][:]
h   = ncv['h'][:]
mcr = ncv['mask_rho'][:]
tnc.close

h    = ma.masked_where(mcr==0,h)
mask = h.mask

h2   = h.copy()
h2[mask] = 0.0

tnci = netCDF4.Dataset(inp_prior)
ncvi = tnci.variables
x    = ncvi ['x_rho'][:]
y    = ncvi ['y_rho'][:]
hp   = ncvi ['h']    [:]
mp   = ncvi ['mask_rho']    [:]

hprior =hp.copy() 

smooth_more = True
every = 30
coef  = 1.0 
print every, coef
if smooth_more:
    hr0 = hp[::every,::every] * mp [::every,::every]
    xr0 = x [::every,::every]
    yr0 = y [::every,::every]
    hp  =interpg(xr0,yr0,hr0,x,y)
    hp  = interp3(xr0,yr0,hr0,x,y,'csa') 
    hp  = hp * coef 

hp_all  = pl.zeros_like(h)
hp_all[j1:j2:k,i1:i2:k] = hp[j1:j2:k,i1:i2:k]
hp_all[mask] = 0.0

alfa = pl.zeros_like(hp_all)
alfa[j1:j2:k,i1:i2:k] = 1.0
#dep_final= h2 * 0.0
nrow = 15

ms1  = linspace(0,1,nrow)

#west
for jm1 in range (j1,j2): 
   alfa[jm1,i1:nrow+i1] = ms1
   alfa[jm1,i1:nrow+i1] = ms1

#east
for im1 in range (j1,j2):  
   alfa[im1,i2-nrow:i2] = linspace(1,0,nrow)

#south
for im1 in range (i1,i2): 
   alfa[j1:nrow+j1,im1] = ms1

#north
for im1 in range (i1,i2): 
   alfa[j2-nrow:j2,im1] = linspace(1,0,nrow)

#imshow(flipud(alfa))
   
dep_f = alfa * hp_all + (1.0-alfa ) * h2


def smoothing(hf):
    h2 = hf.copy()
    print 'smooth the bathymetry and straighten out the edges'
    buf = 4
    for i in range(2, hf.shape[1]-2):
        for j in range(2, hf.shape[0]-2):
            #Saeed invented diagonal smoothing :D
            coef = 1./9.
            h2[j, i] = coef*(hf[j, i]+
                             hf[j+1,i+1]+hf[j+2,i+2]+
                             hf[j+1,i-1]+hf[j+2,i-2]+
                             hf[j-1,i+1]+hf[j-2,i+2]+
                             hf[j-1,i-1]+hf[j-2,i-2])
    return h2

dep0 = dep_f.copy()
buf  = 5
#south
hf   = dep_f[j1-buf:nrow+j1+buf,i1-buf:i2+buf]
dep_f   [j1-buf:nrow+j1+buf,i1-buf:i2+buf] = smoothing(hf)

##north
hf = dep_f[j2-nrow-buf:j2,i1-buf:i2+buf]
dep_f   [j2-nrow-buf:j2,i1-buf:i2+buf]     = smoothing(hf)
#
##east
# j2land=94
# hf=dep_f[j1-buf:j2land+buf,i2-nrow-buf:i2+buf]
# dep_f   [j1-buf:j2land+buf,i2-nrow-buf:i2+buf]=smoothing(hf)
# #
# ##west
# j2land=100
# hf=dep_f[j1-buf:j2land+buf,i1-buf:nrow+i1+buf]
# dep_f   [j1-buf:j2land+buf,i1-buf:nrow+i1+buf]=smoothing(hf)
# #

#dep_f[dep_f < 0.1]= 0.1
#dep_f=dep_f+(mcr-1)*20
#dep_f[dep_f<-5.0]=land


# try:
#     nc1.matlab_inp = read_me 
# except:
#     pass
# nc1.base_dir   = os.getcwd()
# nc1.close()
# 
# comm='cp  '+'tmp.nc  '+new_prior
# os.system(comm)
# 
# if from_mat:
#     nc_name = matfile.replace('.mat','.nc')
#     #to keep track of stuff
#     os.system('cp  '+'tmp.nc  '+nc_name)
#     os.system('mkdir     '+ nc_name[:-3])
#     os.system('mv  *py   '+ nc_name[:-3])
#     os.system('mv  *m    '+ nc_name[:-3])
#     os.system('mv  *mat  '+ nc_name[:-3])
#     os.system('mv  *nc   '+ nc_name[:-3])
#     os.system('mv  *txt  '+ nc_name[:-3])
#     os.system('mv  *inp  '+ nc_name[:-3])
#     os.system('mv  '+ nc_name[:-3]+'/'+new_prior+ '  .')



plt.figure(33)
plt.clf()
ax = plt.subplot(111)
ax.set_aspect(aspect=1)
plt.pcolor(xc,yc,dep_f)
plt.colorbar()

x1 = [xc.min(), xc.min(), xc.max(), xc.max()]
y1 = [yc.max(), yc.min(), yc.min(), yc.max()]
beta1 = [1.0, 1.0, 1.0, 1.0]
grid_tmp2 = BoundaryInteractor(x1,y1,beta1)
plt.draw()

# from pynmd.plotting.points_inside import inside_poly as points_inside_poly 
# verts   = grid_tmp2.verts
# inside  = points_inside_poly(np.vstack( (xc.flat, yc.flat) ).T,verts)
# inside2 = np.reshape(inside,xc.shape)
# jj,ii=dep_f.shape
# deplim = 1.8
# for iy in range(jj):
#     for ix in range (ii):
#         if inside2[iy,ix] and dep_f[iy,ix]< deplim:
#             if not dep_f.mask[iy,ix]:
#                 dep_f[iy,ix] = deplim

dep_f = dep_f * mcr + -10 * np.ones_like(mcr) * (1-mcr)
plt.figure(34)
plt.clf()
ax = plt.subplot(111)
ax.set_title('New Prior')
ax.set_aspect(aspect=1)
plt.pcolor(xc,yc,dep_f)
plt.clim (100,700)
plt.colorbar()
plt.savefig('new_prior.png',dpi=450)

plt.figure(35)
plt.clf()
ax = plt.subplot(111)
ax.set_title('Etopo2 - New Prior')
ax.set_aspect(aspect=1)
plt.pcolor(xc,yc,hprior-dep_f)
plt.clim (-100,100)
plt.colorbar()
plt.savefig('etpo2-new_prior.png',dpi=450)

plt.figure(36)
plt.clf()
ax = plt.subplot(111)
ax.set_title('Etopo2')
ax.set_aspect(aspect=1)
plt.pcolor(xc,yc,hprior)
plt.clim (100,700)
plt.colorbar()
plt.savefig('prior_etopo2.png',dpi=450)


plt.show()

nc1=netCDF4.Dataset('tmp.nc','r+')
ncv1=nc1.variables
ncv1['h'][:]=dep_f[:]
ncv1['h'].missing_value=land
ncv1['h'].valid_min = -1.0
ncv1['h'].valid_max = 5000.0

os.system('cp -vf tmp.nc  ' + new_prior)

