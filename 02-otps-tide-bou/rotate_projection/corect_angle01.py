from numpy import *
import pylab as pl
import netCDF4
import datetime
import os as os
import netcdftime
import sys as sys


print '   >     Rotation according to Geo-projection change will be applied ...'
print '             -->  Different from curvilinear local rotation for ROMS ...'
from mpl_toolkits.basemap import Basemap
proj    = Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)


fngrd0 = 'nri_coastal_tides_2012-05-01_v2.nc'
fngrd  = 'nri_coastal_tides_2012-05-01_v2_rotated.nc'

comm='cp '+fngrd0+'  '+fngrd
os.system(comm)

fgrd = netCDF4.Dataset(fngrd,'r+')
nvargrd=fgrd.variables
latg=nvargrd['lat_rho'][:]
long=nvargrd['lon_rho'][:]
nlat,nlon=long.shape


cang=nvargrd['tide_Cangle'][:]
cangr=pl.zeros_like(cang)


u0=pl.ones_like (long)
v0=pl.zeros_like(latg)

ur=pl.zeros_like (long)
vr=pl.zeros_like(latg)

ur,vr=proj.rotate_vector(u0,v0,long,latg)
alpha=pl.arccos(ur)*180./pi

for id in range(len(cang)):
    cangr[id,:]=cang[id,:]-alpha

cangr[cangr<0.]  =cangr[cangr<0.]  +360.
cangr[cangr>360.]=cangr[cangr>360.]-360.


nvargrd['tide_Cangle'][:]=cangr

#dimensions:
#    two = 2 ;
#    eta_rho = 300 ;
#    xi_rho = 240 ;
#    tide_period = UNLIMITED ; // (8 currently)
#double tide_Cangle(tide_period, eta_rho, xi_rho) ;
#        tide_Cangle:long_name = "tidal current inclination angle" ;
#        tide_Cangle:units = "degrees between semi-major axis and East" ;
#        tide_Cangle:field = "tide_Cangle, scalar" ;


angn = fgrd.createVariable('rotation_ang', 'float', ('eta_rho', 'xi_rho',))
angn.long_name = 'rotation of grid base on projection change'
angn.units = 'deg'
angn.valid_min = 0.0
angn.valid_max = 360
angn[:] = alpha


fgrd.close()

