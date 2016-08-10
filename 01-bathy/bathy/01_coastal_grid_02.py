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

from octant.grid import *
#from grid import *
import netCDF4
import octant.csa as csa
import roms as roms
import cPickle
import datetime,sys
from okean import calc
import numpy as np
import matplotlib.pyplot as plt
import string
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap

inp_dir = 'inp/'
vv = True

smith_9_1_bathy  = inp_dir + 'data_v9-1_smith_nri.nc'
gebco_bathy      = inp_dir + 'GEBCO_2014_2D_nri.nc'
etopo2_bathy     = inp_dir + 'ETOPO2v2c_f4_nri.nc'

topo  = 'gebco'
#topo    = 'smith'
#topo   = 'etopo2'

print topo
if   topo == 'gebco':
    bathy_inp = gebco_bathy
    vvv       = 'gebco_'
    llon      = 'lon'
    llat      = 'lat'
    lbat      = 'elevation'
elif topo == 'smith':
    bathy_inp = smith_9_1_bathy 
    vvv       = 'smith_9_1_'
    llon      = 'longitude'
    llat      = 'latitude'
    lbat      = 'topo'
elif topo == 'etopo2':
    bathy_inp = etopo2_bathy 
    vvv       = 'etopo2_'
    llon      = 'x'
    llat      = 'y'
    lbat      = 'z'
else:
    print '  >>>  choose correct option >>'

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

#read GEBCO
ncl  = netCDF4.Dataset(bathy_inp, 'r')
ncv  = ncl.variables
lon_geb  =  ncv[llon][:]
lat_geb  =  ncv[llat][:]
bat_geb  = -ncv[lbat][:]
bat_geb  =  np.ma.masked_where(bat_geb < -2,bat_geb )  

lon_geb2,lat_geb2 = np.meshgrid(lon_geb,lat_geb)

# extent of new grid
lonlim1 = -79.0 ,-75.0
latlim1 =  31.0 , 35.1
limits = np.array([lonlim1,latlim1]).flatten()
    
# simplify data:
print 'Simplify data ...'
i0,i1,j0,j1 = calc.ij_limits(lon_geb2,lat_geb2,lonlim1,latlim1,margin=1)
lon_geb = lon_geb2 [j0:j1,i0:i1]
lat_geb = lat_geb2 [j0:j1,i0:i1]
bat_geb = bat_geb [j0:j1,i0:i1]

#read altimeter data positions
lon_sat,lat_sat = read_path(limits)

#from below plots
lon_tmp = np.array([-78.834167780118207, -77.217793492727196, -75.11248036395456, -76.700093543029013])
lat_tmp = np.array([ 33.94159435781355  , 31.065483526156946 , 32.376990065392356, 35.138056463782696])
#GRID    = np.array([lon_tmp,lat_tmp, [1.0, 1.0, 1.0, 1.0]])

if vv:
    x1 = [lonlim1[0], lonlim1[0], lonlim1[1], lonlim1[1]]
    y1 = [latlim1[1], latlim1[0], latlim1[0], latlim1[1]]
    x1 = lon_tmp
    y1 = lat_tmp
    beta1 = [1.0, 1.0, 1.0, 1.0]
    plt.figure(1)
    plt.clf()
    ax = plt.subplot(111)
    ax.set_aspect(aspect=1)
    plt.contourf(lon_geb,lat_geb,bat_geb,N=50)
    plt.scatter(lon_sat,lat_sat,s=10,c='k')
    grid_tmp2 = BoundaryInteractor(x1,y1,beta1)
    plt.draw()
    plt.savefig('pic/fig01_'+vvv+'.png',dpi=450)

#tmerc            Transverse Mercator                     
#width            width of desired map domain in projection coordinates (meters).
#height           height of desired map domain in projection coordinates (meters).
#lon_0            center of desired map domain (in degrees).
#lat_0            center of desired map domain (in degrees).
#rsphere          radius of the sphere used to define map projection

# to UTM
#proj = Basemap(llcrnrlon=-77.5,llcrnrlat =34.3,
#                  urcrnrlon=-77.1 ,urcrnrlat=34.9,
#                  resolution='i',projection='tmerc',
#                  lon_0=-77.33,lat_0=34.53)

if False:
    proj = pyproj.Proj(proj='utm',zone='18',ellps='WGS84')
    xu_orig , yu_orig  = proj(-77.3382, 34.5279)
    xg = np.array([-147439.98467231489, 2892.5195482667768, 214783.59230744862, 59721.210069045541])  + xu_orig
    yg = np.array([-46597.880410829792,-384198.02997515333,-263672.83262589388, 66502.105108662974])  + yu_orig
else:
    proj    = Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)
    xu_orig , yu_orig  = proj(-77.3382, 34.5279)
    xg = np.array([-191112.42519508279, -248001.34775263432, 69124.1354405249  , 105436.2136687492])  + xu_orig
    yg = np.array([47756.549282671302 , -363780.33730387117,-392829.99988645065, 23548.49713052169])  + yu_orig

x_geb   , y_geb    = proj(lon_geb, lat_geb)
x_sat   , y_sat    = proj(lon_sat, lat_sat)
#xg      ,  yg      = proj(lon_tmp,lat_tmp)

#
x_geb_o = x_geb - xu_orig
y_geb_o = y_geb - yu_orig

x_sat_o = x_sat - xu_orig
y_sat_o = y_sat - yu_orig

xg_o      = xg - xu_orig
yg_o      = yg - yu_orig

if vv==True:
    plt.figure(2)
    plt.clf()
    ax = plt.subplot(111)
    ax.set_aspect(aspect=1)
    plt.contourf(x_geb_o,y_geb_o,bat_geb,N=50)
    plt.scatter(0, 0          , s=20, c='k')
    plt.scatter(x_sat_o,y_sat_o, s=10, c='r')
    plt.plot(xg_o , yg_o ,'-k')
    grid_tmp2 = BoundaryInteractor(xg_o,yg_o,beta1)
    plt.draw()

long,latg = proj(xg,yg,inverse=True)
GRID      = np.array([long,latg, [1.0, 1.0, 1.0, 1.0]])

print 'Gridgen  ...'
# focus the grid near the estuary mouth (parameters were determined
# by trial-and-error)
#if False:
#    def focus(x, y, xo=1.5, yo=1.5):
#        xf = np.tan((x - xo)*2.0)
#        yf = np.tan((y - yo)*2.0)
#        xf -= xf.min()
#        xf /= xf.max()
#        yf -= yf.min()
#        yf /= yf.max()
#        return xf, yf

ul_idx = 0  # index of upper left hand corner
#shp = (201,201)        # shape of the large grid
shp = (191,151)        # shape of the large grid

zooming = False
if zooming:
        foc = octant.grid.Focus()
        foc.add_focus_x(0.45  ,  factor=4 , Rx=0.5)
        foc.add_focus_y(0.5   ,  factor=4 , Ry=0.5)
        grdf=Gridgen(GRID[0],GRID[1], GRID[2], shp,proj=proj, focus=foc, ul_idx=ul_idx) #, verbose=True)
    
else:
        grdf=Gridgen(GRID[0],GRID[1], GRID[2], shp,proj=proj, focus=None, ul_idx=ul_idx) #, verbose=True)

####>>>>>#### copy and paste from grid.py 
grdf.lon_rho, grdf.lat_rho = grdf.proj(grdf.x_rho , grdf.y_rho , inverse=True)
grdf.lon_u,   grdf.lat_u =   grdf.proj(grdf.x_u   , grdf.y_u   , inverse=True)
grdf.lon_v,   grdf.lat_v =   grdf.proj(grdf.x_v   , grdf.y_v   , inverse=True)
grdf.lon_psi, grdf.lat_psi = grdf.proj(grdf.x_psi , grdf.y_psi , inverse=True)
grdf.lon_vert,grdf.lat_vert= grdf.proj(grdf.x_vert, grdf.y_vert, inverse=True)

grdf.x_rho = grdf.x_rho - xu_orig
grdf.y_rho = grdf.y_rho - yu_orig
 
grdf.x_u   = grdf.x_u   - xu_orig
grdf.y_u   = grdf.y_u   - yu_orig
 
grdf.x_v   = grdf.x_v   - xu_orig
grdf.y_v   = grdf.y_v   - yu_orig
 
grdf.x_psi = grdf.x_psi - xu_orig
grdf.y_psi = grdf.y_psi - yu_orig
 
grdf.x_vert= grdf.x_vert- xu_orig
grdf.y_vert= grdf.y_vert- yu_orig

grdf.f = 2.0 * 7.29e-5 * np.cos(grdf.lat_rho * np.pi / 180.0)
######<<<<####

dx1 = abs(grdf.x_rho[:-1,:-1])-abs(grdf.x_rho[1:,1:])
dx1 = abs(dx1)
dy1 = grdf.y_rho[:-1,:-1]-grdf.y_rho[1:,1:]
dd1 = np.sqrt(dx1**2+dy1**2)


print 'grid shape >  ', grdf.x_rho.shape
print 'min and max grid size (dd) >  ', dd1.min(),dd1.max()
print 'min and max grid size (dx) >  ', dx1.min(),dx1.max()
print 'min and max grid size (dy) >  ', dy1.min(),dy1.max()


if vv:
    plt.figure(22)
    plt.clf()
    ax = plt.subplot(111)
    ax.set_aspect(aspect=1)
    img4 = ax.contourf(x_geb_o,y_geb_o,bat_geb,N=50)
    ax.plot(grdf.x_rho  ,grdf.y_rho  ,'-k',alpha=0.2)
    ax.plot(grdf.x_rho.T,grdf.y_rho.T,'-k',alpha=0.2)
    ax.scatter(x_sat_o,y_sat_o, s=10, c='r')
    ax.plot(xg-xu_orig,yg-yu_orig,'-k')
    cb=plt.colorbar(img4)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    ax.set_xlabel('Local X (m)')
    ax.set_ylabel('Local Y (m)')
    plt.savefig('pic/fig022'+vvv+'.png',dpi=450)
    #plt.close()  

if vv:
    vmax = 3000
    levels = np.arange(0,vmax,vmax/50) 
    
    
    plt.figure(23)
    plt.clf()
    ax = plt.subplot(111)
    ax     = plt.subplot(111,aspect=(1.0/np.cos(np.mean(grdf.lat_rho)*np.pi/180.0)))
    img4 = plt.contourf(lon_geb, lat_geb,bat_geb,levels=levels ,shading='faceted',cmap=plt.cm.jet,extend='both')
    ax.plot(grdf.lon_rho  ,grdf.lat_rho  ,'-k',alpha=0.2)
    ax.plot(grdf.lon_rho.T,grdf.lat_rho.T,'-k',alpha=0.2)
    ax.scatter(lon_sat, lat_sat, s=10, c='r')
    cb=plt.colorbar(img4)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    plt.gca().patch.set_facecolor('0.5')

    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    plt.savefig('pic/fig022_lon_lat'+vvv+'.png',dpi=450)
    #plt.close()  









#sys.exit()

xf = x_geb_o[~bat_geb.mask]
yf = y_geb_o[~bat_geb.mask]
bf = bat_geb[~bat_geb.mask]


# Interpolate the bathymetry points to the curvelinear grid
# The bathymetry points:
print 'interp ...'

interpm='csa'
if interpm=='tri':
    print 'tri interp ...'
    from delaunay import  *
    tri=triangulate.Triangulation(xf, yf)
    interp_b = tri.nn_extrapolator(bf)
    b_gridf  = interp_b(grdf.x_rho.flatten(),grdf.y_rho.flatten())
    b_grid   = b_gridf.reshape(grdf.x_rho.shape)
elif interpm == 'csa' :
    print 'csa interp ...'
    import octant.csa as csa
    csa_interp = csa.CSA(xf, yf,bf)
    b_grid = csa_interp(grdf.x_rho,grdf.y_rho)
elif interpm == 'grid' :
    print 'griddata interp ...'
    b_grid = pl.griddata(xf,yf,bf,grdf.x_rho,grdf.y_rho)

hf = b_grid.copy()

print 'modifying min/max depth'
h_min =  5.0
h_max =  10000.0

hf[hf< 0 ]              = np.nan
hf[ (hf>0) &(hf<h_min)] = h_min
hf[hf>h_max]            = h_max

hf = np.ma.masked_where(np.isnan(hf),hf)

if True :
   fig=plt.figure()
   ax = plt.subplot(111)
   ax.set_aspect(aspect=1)
   img3=plt.pcolor(grdf.x_rho,grdf.y_rho,hf)
   cb=plt.colorbar(img3)
   plt.clim (0,10)
   grid_tmp2 = BoundaryInteractor(xg_o,yg_o,beta1)
   plt.draw()

#sys.exit() 
verts1 = [(-191112.42519508279, 47756.549282671418),
 (-209572.74899396364, 11745.472837022215),
 (-192595.88782696164, 8350.1006036218023),
 (-164125.14801372297, 3431.057595508435),
 (-132060.77138476784, -17428.0261254785),
 (-118442.03077354498, -28633.319033446678),
 (-109994.96381215358, -27771.373425141428),
 (-97238.168809235969, -27254.206060158293),
 (-77068.641574893234, -9842.9047723923431),
 (-52253.835513078375, 3822.9376257546246),
 (-31040.746091393172, 5327.3379337799706),
 (44810.467439468368, 1879.5555005589849),
 (87390.580489747459, -22254.921531987871),
 (108077.27508907334, -29495.264641751928),
 (118646.5669014086, 30985.915492957865)]

inside  = points_inside_poly(np.vstack( (grdf.x_rho.flatten(),grdf.y_rho.flatten()) ).T,verts1)
inside2 = np.reshape(inside,grdf.x_rho.shape)
hf[inside2] = np.nan

# ################# > take thi one out >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# verts1 = [(-174487.23591549281, 16838.531187122804),
#  (-214099.91197183088, -74270.623742454685),
#  (-207875.06287726347, 37210.76458752522),
#  (122041.93913480907, 41172.032193159044),
#  (98840.228873239597, -146705.23138832994),
#  (95444.856639839243, -7494.9698189133778)]
# 
# inside  = points_inside_poly(np.vstack( (grdf.x_rho.flatten(),grdf.y_rho.flatten()) ).T,verts1)
# inside2 = np.reshape(inside,grdf.x_rho.shape)
# hf[inside2] = np.nan

###  Create mask for grid
hf = np.ma.masked_where(np.isnan(hf),hf)
grdf.mask[hf.mask]  = 0
grdf.mask[~hf.mask] = 1

hf[np.isnan(hf)] = -2.0
grdf.h = hf
#############################################
if vv :
    fig=plt.figure()
    fig.clf()
    ax = plt.subplot(111)
    ax.set_aspect(aspect=1)
    img3=plt.pcolor(grdf.lon_rho,grdf.lat_rho,grdf.mask_rho)
    #gca().xlabel('Lon')
    #gca().ylabel('Lat')
    plt.savefig('pic/mask'+vvv+'.png',dpi=450)
    plt.close()    
    ##
    fig  = plt.figure()
    fig.clf()
    vmax = grdf.h.max()
    vmax = 2500
    h_masked = np.ma.masked_where(grdf.mask_rho==0,grdf.h )
    levels = np.arange(0,vmax,vmax/50)   
    ax     = plt.subplot(111,aspect=(1.0/np.cos(np.mean(grdf.lat_rho)*np.pi/180.0)))
    img3 = plt.contourf(grdf.lon_rho,grdf.lat_rho,h_masked , levels = levels ,shading='faceted',cmap=plt.cm.jet,extend='both')
    cb   = plt.colorbar(img3)
    color_lab = 'Depth (m)'
    cb.set_label(color_lab)
    plt.gca().patch.set_facecolor('0.5')
    plt.xlabel('Lon')
    plt.ylabel('Lat')
    plt.savefig('pic/h'+vvv+'.png',dpi=450)







if True:
    print 'Writing netcdf file..'
    rom_grd_name='nri_coastal_grd_roms_v2'+vvv+'.nc'
    roms.write_grd(grdf, rom_grd_name, verbose=True)
    
"""
# saeed moghimi
# moghimis@gmail.com
# osu

vv=False


inc_local=True
#read local bathymetry file
local_raw=False

if inc_local:
    if local_raw:
        print 'Read local 2012 bathymetry raw data ..'
        # saeed moghimi
        # moghimis@gmail.com
        #2012-07-16 21:17:03 
        
        fname='LARCsurvey_1_20120416XXXXXXXXXXXXXXXXXXXXXXXXXXX.txt'
        fin1 = open(fname,'r')
        count_lines = len(fin1.readlines())
        fin1.close()
        
        #restart reading the file
        fin1 = open(fname,'r')
        
        latl = zeros([count_lines], float)
        lonl = zeros([count_lines], float)
        elhl = zeros([count_lines], float)
        xl  = zeros([count_lines], float)
        yl  = zeros([count_lines], float)
        bl  = zeros([count_lines], float)
        
        
        for i in range(count_lines):
           a2= fin1.readline()
           latl[i],lonl[i],elhl[i],xl[i],yl[i],bl[i]= fromstring(a2, dtype=float, sep=',')        #Lat,Lon,Ellipsoid HT,XS,AS,NAVD88 HT
           
              
        fin1.close()
    
    else:
        print 'Read local 2012 bathymetry falk gridded data ...'
        #ncl= netCDF4.Dataset('LARCsurvey_1_gridded_asxsz_20120416.nc')
        ncl= netCDF4.Dataset('LARCsurvey_3_gridded_asxsz_20120510.nc')
        nvarl=ncl.variables
        lonl  = nvarl['lonc'][:]
        latl  = nvarl['latc'][:]
        bl = nvarl['bathymetry'][:]
        maskl=bl.mask

    
    x1 = [lonl.min(), lonl.min(), lonl.max(), lonl.max()]
    y1 = [latl.max(), latl.min(), latl.min(), latl.max()]
    beta1 = [1.0, 1.0, 1.0, 1.0]
    levels=linspace(0,20,41)
    if vv==True:
        pl.figure(23)
        pl.clf()
        ax = pl.subplot(111)
        ax.set_aspect(aspect=1)
        if local_raw: pl.plot(lonl,latl,'.')
        if not local_raw:  pl.pcolor(lonl,latl,bl)
        grid_tmp2 = BoundaryInteractor(x1,y1,beta1)
        pl.draw()
        pl.close()

    # convert from NAVD to MSL (from Gerg)
    navd2msl=-1.588+1.453;
    bl_temp=bl * -1.0
    bl=bl_temp+navd2msl

#Cut1
#vertsl=[(-77.341831760190857, 34.539059671945694),
# (-77.345285758738726, 34.535595305429858),
# (-77.346001452131532, 34.533863122171944),
# (-77.345752515299253, 34.532378393665155),
# (-77.343823254849084, 34.528277714932123),
# (-77.344881236386271, 34.526404128959271),
# (-77.349175396743078, 34.523540723981895),
# (-77.348770874390624, 34.522515554298636),
# (-77.343543200912777, 34.516223133484161),
# (-77.33377243024583, 34.518556278280542),
# (-77.331158593506899, 34.518909785067869),
# (-77.3271444870864, 34.521631787330314),
# (-77.324592884555543, 34.524000282805424),
# (-77.32129447152785, 34.526722285067869),
# (-77.323534903018356, 34.529408936651578),
# (-77.327953531791309, 34.534782239819002),
# (-77.332870034228819, 34.53085831447963),
# (-77.334737060470914, 34.532024886877821),
# (-77.337599834042109, 34.532343042986419),
# (-77.341084949694022, 34.533686368778277),
# (-77.342765273311898, 34.534570135746598),
# (-77.341396120734373, 34.535312499999996),
# (-77.339809148428586, 34.536408371040721)]



#Cut smaller
vertsl=[(-77.341905051816951, 34.53831649001765),
 (-77.343609657805416, 34.536514423076916),
 (-77.344846931561079, 34.535065045248864),
 (-77.345377191742074, 34.533898472850673),
 (-77.345518594457005, 34.533120757918546),
 (-77.345377191742074, 34.532449095022621),
 (-77.345214749574822, 34.531839095867191),
 (-77.344811580882336, 34.530893665158366),
 (-77.344508435815641, 34.53037938076487),
 (-77.343550988275396, 34.528637140158871),
 (-77.343256151018082, 34.527606052036198),
 (-77.343680359162889, 34.526616233031667),
 (-77.344630931016098, 34.525911528848063),
 (-77.346810056424317, 34.524591776840275),
 (-77.348145154385691, 34.523763095347014),
 (-77.346932621606328, 34.521808540723974),
 (-77.343220800339353, 34.517566459276011),
 (-77.341806773190029, 34.517425056561081),
 (-77.339968537895913, 34.51784926470588),
 (-77.33413567590496, 34.519616798642531),
 (-77.326393877262433, 34.523293269230763),
 (-77.325859173038012, 34.524571753055746),
 (-77.324838447398179, 34.525166855203615),
 (-77.322823458710403, 34.526934389140266),
 (-77.32819676187782, 34.533191459276011),
 (-77.332332791289588, 34.530186651583705),
 (-77.333676117081438, 34.530115950226239),
 (-77.335372949660623, 34.531141119909499),
 (-77.33820100395927, 34.531529977375563),
 (-77.340498798076908, 34.532519796380086),
 (-77.342825394142636, 34.533733185235704),
 (-77.343377599538059, 34.534800782333505),
 (-77.34219956136117, 34.535739531505712),
 (-77.341444880654109, 34.536052447896445),
 (-77.340137994551625, 34.536383771133693),
 (-77.341187184802905, 34.537138451840761),
 (-77.341518508040153, 34.537948353087373)]

####################################################33

print 'Read regional 2011 netcdf  ...'

#ncks -O -d lonc,1550,2300  -d latc,900,1600  topo10m_nri.nc topo10m_nri_small.nc
nc1= netCDF4.Dataset('topo10m_nri_small.nc')
#nc1= netCDF4.Dataset('topo10m_nri_lat_50m.nc')
#nc1= netCDF4.Dataset('topo10m_nri_lat_100m.nc')

#for final
#('topo10m_nri_small.nc') 

# for prior('topo10m_nri_lat_50m.nc')
#            topo10m_nri_lat_100m.nc 
nvar=nc1.variables


ist=0
ien=-1
inn=1

jst=0
jen=-1
jnn=1


lon = nvar['lonc'][ist:ien:jnn]
lat = nvar['latc'][jst:jen:inn]
b = nvar['bathymetry']    [jst:jen:inn,ist:ien:jnn]
b= ma.masked_where(b<=-30.0,b)

lon,lat=meshgrid(lon,lat)

# extent of new grid
lonlim1=-77.37,-77.30
latlim1=34.50,34.55

# simplify data:
if True:    
    print 'Simplify data ...'
    i0,i1,j0,j1=calc.ij_limits(lon,lat,lonlim1,latlim1,margin=3)
    lon=lon[j0:j1,i0:i1]
    lat=lat[j0:j1,i0:i1]
    b=b[j0:j1,i0:i1]

subt=False
if inc_local and True:
    subt=True
    print 'Subtract local from regional bathy  ...'
    inside = points_inside_poly(np.vstack( (lon.flatten(), lat.flatten()) ).T,vertsl)
    inside=pl.reshape(inside,lon.shape)
    b[inside]=0.000911
    b=ma.masked_where(b==0.000911,b)
    lon=ma.masked_where(b==0.000911,lon)
    lat=ma.masked_where(b==0.000911,lat)
    maskb=b.mask


#####################################
lonmin=lon.min()-0.1
lonmax=lon.max()+0.1
lonmean= (lonmin + lonmax)/2

latmin=lat.min()-0.1
latmax=lat.max()+0.1
latmean= (latmin + latmax)/2

#   llcrnrlon        longitude of lower left hand corner of the desired map domain (degrees).
#   llcrnrlat        latitude of lower left hand corner of the desired map  domain (degrees).
#   urcrnrlon        longitude of upper right hand corner of the desired map domain (degrees).
#   urcrnrlat        latitude of upper right hand corner of the desired map  domain (degrees).

#m=Basemap(projection='merc',lon_0=13.5,lat_ts=54.2,resolution='i')
#boundary resolution must be one of 'c','l','i','h' or 'f'

#proj = Basemap(projection='merc',
#    resolution=None,
#    llcrnrlon=lonmin ,
#    llcrnrlat=latmin ,
#    urcrnrlon=lonmax ,
#    urcrnrlat=latmax ,
#    lon_0=lonmean,
#    lat_0=latmean)

proj    = Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)

proj_str="Basemap(resolution=None,projection='npstere',lon_0=-45,boundinglat=70)"

#proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
#####################################


x1 = [lon.min(), lon.min(), lon.max(), lon.max()]
y1 = [lat.max(), lat.min(), lat.min(), lat.max()]
beta1 = [1.0, 1.0, 1.0, 1.0]

levels=linspace(0,15,61)

if vv:
    pl.figure(2)
    pl.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    pl.contourf(lon,lat,b,levels)
    pl.colorbar()
    grid_tmp2 = BoundaryInteractor(x1,y1,beta1)
    pl.draw()


add_ch=True

if add_ch:
    verts=[(-77.356835434254165, 34.551136047659341),
        (-77.357353864886065, 34.550594730196956),
        (-77.357262377127498, 34.550354144658122),
        (-77.358604197586516, 34.549030924194525),
        (-77.358970148620784, 34.548850485040397),
        (-77.359183620057451, 34.547316752230323),
        (-77.358878660862217, 34.546655141998521),
        (-77.358452698354156, 34.545759956283177),
        (-77.3588481649427, 34.54512140918844),
        (-77.359305603735535, 34.543106505300692),
        (-77.359125899598865, 34.542495746103242),
        (-77.359366595574585, 34.537422671945698),
        (-77.359671554769818, 34.537091866829797),
        (-77.360189985401718, 34.535197255711466),
        (-77.361409822182637, 34.533092132246651),
        (-77.361653789538821, 34.53215986328366),
        (-77.36268319348514, 34.530723891659811),
        (-77.364635753393458, 34.528223751728525),
        (-77.36558237369664, 34.527140374087054),
        (-77.366322592834152, 34.525130938052456),
        (-77.367753683166669, 34.522875448625868),
        (-77.368920901387028, 34.521189882646745),
        (-77.369579557039202, 34.520209870212632),
        (-77.370680872884378, 34.51875295162629),
        (-77.372589781531772, 34.517831354090049),
        (-77.374859786886816, 34.516929158319414),
        (-77.390503084659628, 34.511885063783595),
        (-77.390256344947119, 34.511392956999615),
        (-77.374613047174307, 34.516437051535434),
        (-77.372096302106769, 34.517462274002064),
        (-77.37007303646422, 34.518364469772699),
        (-77.368789989959183, 34.520004825719305),
        (-77.368148466706671, 34.520989039287272),
        (-77.367112159914157, 34.52246535963922),
        (-77.36568106958164, 34.524802866863133),
        (-77.364916852927792, 34.526746688659863),
        (-77.364100109894267, 34.528500675708123),
        (-77.362385691607372, 34.530836642820063),
        (-77.361501309941204, 34.532039570514243),
        (-77.360891391550751, 34.532610961168977),
        (-77.359549571091733, 34.534956670172626),
        (-77.358970148620784, 34.536971574060381),
        (-77.359214115976982, 34.537693330676888),
        (-77.358970148620784, 34.542625334223018),
        (-77.358665189425551, 34.543016285723631),
        (-77.358055271035099, 34.54485075045725),
        (-77.358360230230332, 34.54596345857437),
        (-77.358695685345083, 34.546745361575589),
        (-77.358360230230332, 34.548369313962731),
        (-77.35839072614985, 34.549121143771586),
        (-77.357079401610349, 34.550293998273418),
        (-77.356743946495598, 34.55020377869635),
        (-77.35631700362228, 34.550594730196956),
        (-77.355951052587997, 34.550865388928152),
        (-77.356499979139414, 34.551436779582886)]
    
    inside = points_inside_poly(np.vstack( (lon.flatten(), lat.flatten()) ).T,verts)
    inside2=pl.reshape(inside,lon.shape)
    b[inside2]=3.0
    
    verts=[(-77.397657479560181, 34.510317234053581),
     (-77.397073816226609, 34.509245919100259),
     (-77.379778950079526, 34.514742994188602),
     (-77.379901826570801, 34.515252307854929)]
    inside = points_inside_poly(np.vstack( (lon.flatten(), lat.flatten()) ).T,verts)
    inside2=pl.reshape(inside,lon.shape)
    b[inside2]=3.0
    
    ### add more prominent channel to snear_frey area to have real channel in numerical  
    verts=[(-77.402838464911852, 34.580271672788847),
     (-77.401119457237087, 34.577154280113817),
     (-77.401641758078412, 34.574432382813072),
     (-77.397968940425301, 34.575966115623153),
     (-77.39368573729746, 34.578233974794941),
     (-77.395492196695031, 34.580821851247208),
     (-77.395767466698473, 34.582635402461783),
     (-77.399019093614115, 34.581310898765743)]
    inside = points_inside_poly(np.vstack( (lon.flatten(), lat.flatten()) ).T,verts)
    inside2=pl.reshape(inside,lon.shape)
    b[inside2]=3.4


#
xu, yu = proj(lon, lat)
xu_orig,yu_orig=proj(-77.3382, 34.5279)

#grid 1
long=(-77.262317493920179,
  -77.608511480370481,
  -77.396610679993216,
  -77.035920840024531)

latg=(34.847728163203989,
      34.665728368207496,
      34.184204222897264,
      34.395940436800927)

xg,yg=proj(long,latg)
GRID =array([xg,yg, [1.0, 1.0, 1.0, 1.0]])


####################
#grid 2

#xg=array([-1.326e6, -1.343e6, -1.343e6, -1.326e6])
#yg=array([-3.383e6, -3.383e6, -3.423e6, -3.423e6])

#xg=array([-1.314e6, -1.355e6, -1.355e6, -1.314e6])
#yg=array([-3.38e6, -3.38e6, -3.449e6, -3.449e6])


xg=array([-1335000., -1339000., -1339000., -1335000.])
yg=array([-3412000., -3412000., -3416000., -3416000.])


xlen=np.abs(xg[1]-xg[0])/1000.0
ylen=np.abs(yg[2]-yg[0])/1000.0


print 'xlen= ', xlen,'ylen= ',ylen,'  (km)'

#for plotting
GRID=array([xg,yg, [1.0, 1.0, 1.0, 1.0]])

#[grid_tmp.x,grid_tmp.y,grid_tmp.beta]

if vv:
    pl.figure(3)
    pl.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img=pl.contourf(xu,yu,b,levels)
    cb=pl.colorbar(img)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    #pl.gca().set_xlabel('Lon. (Deg.)')
    #pl.gca().set_ylabel('Lat. (Deg.)')
    pl.gca().set_xlabel('Local X (m)')
    pl.gca().set_ylabel('Local Y (m)')
    grid_tmp = BoundaryInteractor(GRID[0],GRID[1], GRID[2])
    #pl.gca().set_ylim( 34.50, 34.55)
    #pl.gca().set_xlim(-77.37,-77.30)
    pl.draw()
    pl.savefig('pic/fig3.png',dpi=450)
    #pl.close()

#for gridgen
long,latg=proj(xg,yg,inverse=True)
GRID=array([long,latg, [1.0, 1.0, 1.0, 1.0]])

#sys.exit()

print 'Gridgen  ...'
# focus the grid near the estuary mouth (parameters were determined
# by trial-and-error)
#if False:
#    def focus(x, y, xo=1.5, yo=1.5):
#        xf = np.tan((x - xo)*2.0)
#        yf = np.tan((y - yo)*2.0)
#        xf -= xf.min()
#        xf /= xf.max()
#        yf -= yf.min()
#        yf /= yf.max()
#        return xf, yf

ul_idx = 1  # index of upper left hand corner
shp = (201,201)        # shape of the large grid
#shp = (101,101)        # shape of the large grid

zooming=False
if zooming:
        foc = octant.grid.Focus()
        foc.add_focus_x(0.45  ,  factor=4 , Rx=0.5)
        foc.add_focus_y(0.5   ,  factor=4 , Ry=0.5)
        grdf=Gridgen(GRID[0],GRID[1], GRID[2], shp,proj=proj, focus=foc, ul_idx=ul_idx) #, verbose=True)
    
else:
        grdf=Gridgen(GRID[0],GRID[1], GRID[2], shp,proj=proj, focus=None, ul_idx=ul_idx) #, verbose=True)

####>>>>>#### copy and paste from grid.py 
grdf.lon_rho, grdf.lat_rho = grdf.proj(grdf.x_rho , grdf.y_rho , inverse=True)
grdf.lon_u,   grdf.lat_u =   grdf.proj(grdf.x_u   , grdf.y_u   , inverse=True)
grdf.lon_v,   grdf.lat_v =   grdf.proj(grdf.x_v   , grdf.y_v   , inverse=True)
grdf.lon_psi, grdf.lat_psi = grdf.proj(grdf.x_psi , grdf.y_psi , inverse=True)
grdf.lon_vert,grdf.lat_vert= grdf.proj(grdf.x_vert, grdf.y_vert, inverse=True)


xu = xu - xu_orig
yu = yu - yu_orig

grdf.x_rho = grdf.x_rho - xu_orig
grdf.y_rho = grdf.y_rho - yu_orig

grdf.x_u   = grdf.x_u   - xu_orig
grdf.y_u   = grdf.y_u   - yu_orig

grdf.x_v   = grdf.x_v   - xu_orig
grdf.y_v   = grdf.y_v   - yu_orig

grdf.x_psi = grdf.x_psi - xu_orig
grdf.y_psi = grdf.y_psi - yu_orig

grdf.x_vert= grdf.x_vert- xu_orig
grdf.y_vert= grdf.y_vert- yu_orig


grdf.f = 2.0 * 7.29e-5 * np.cos(grdf.lat_rho * np.pi / 180.0)
######<<<<####

dx1=abs(grdf.x_rho[:-1,:-1])-abs(grdf.x_rho[1:,1:])
dx1=abs(dx1)
dy1=grdf.y_rho[:-1,:-1]-grdf.y_rho[1:,1:]
dd1=sqrt(dx1**2+dy1**2)


print 'grid shape >  ', grdf.x_rho.shape
print 'min and max grid size (dd) >  ', dd1.min(),dd1.max()
print 'min and max grid size (dx) >  ', dx1.min(),dx1.max()
print 'min and max grid size (dy) >  ', dy1.min(),dy1.max()

#pl.figure()
#pl.pcolor(dd1)
#pl.colorbar()
#pl.clim(0,50)
#pl.show()

if vv:
    pl.figure(4)
    pl.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.contourf(lon,lat,b,levels)
    pl.plot(grdf.lon_rho,grdf.lat_rho,'-k',alpha=0.3)
    pl.plot(grdf.lon_rho.T,grdf.lat_rho.T,'-k',alpha=0.3)
    pl.draw()
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('Lon. (Deg.)')
    pl.gca().set_ylabel('Lat. (Deg.)')
    #pl.gca().set_ylim( 34.4, 34.6)
    #pl.gca().set_xlim(-77.5,-77.2)
    #pl.gca().set_xlabel('Local X (m)')
    #pl.gca().set_ylabel('Local Y (m)')
    mask1=b.mask
    pl.savefig('pic/fig4.png',dpi=450)
    #pl.close()    


if vv:
    pl.figure(444)
    pl.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.contourf(xu,yu,b,levels)
    pl.plot(grdf.x_rho  ,grdf.y_rho  ,'-k',alpha=0.2)
    pl.plot(grdf.x_rho.T,grdf.y_rho.T,'-k',alpha=0.2)
    pl.draw()
    pl.gca().set_ylim( yu.min(), yu.max())
    pl.gca().set_xlim( xu.min(), xu.max())
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('Local X (m)')
    pl.gca().set_ylabel('Local Y (m)')
    mask1=b.mask
    pl.savefig('pic/fig44.png',dpi=450)
    #pl.close()  

#sys.exit()
#
### Interpolation

if subt:
    xf=xu.data[~maskb].flatten()
    yf=yu.data[~maskb].flatten()
    bf= b.data[~maskb].flatten()
else:
    xf=xu.flatten()
    yf=yu.flatten()
    bf= b.flatten()

if inc_local:
    x_local,y_local=grdf.proj(lonl, latl)
    
    x_local=x_local-xu_orig
    y_local=y_local-yu_orig
    
    
    if local_raw:
        lonlmf=x_local
        latlmf=y_local
        blmf=     bl
    else:
        lonlmf=x_local[~maskl].flatten()
        latlmf=y_local[~maskl].flatten()
        blmf=     bl[~maskl].flatten()
    
    xf=np.hstack((xf,lonlmf))
    yf=np.hstack((yf,latlmf))
    bf=np.hstack((bf,blmf))

#print 'tri'
#tri=triangulate.Triangulation(xf, yf)
#print 'interp'
#interp_b=tri.nn_extrapolator(bf)


# Interpolate the bathymetry points to the curvelinear grid
# The bathymetry points:
print 'interp ...'

interpm=csa
if interpm==tri:
    print 'tri interp ...'
    from delaunay import  *
    tri=triangulate.Triangulation(xf, yf)
    interp_b=tri.nn_extrapolator(bf)
    b_gridf = interp_b(grdf.x_rho.flatten(),grdf.y_rho.flatten())
    b_grid=b_gridf.reshape(grdf.x_rho.shape)
elif interpm==csa :
    print 'csa interp ...'
    import octant.csa as csa
    csa_interp = csa.CSA(xf, yf,bf)
    b_grid = csa_interp(grdf.x_rho,grdf.y_rho)
elif interpm==grid :
    print 'griddata interp ...'
    b_grid=pl.griddata(xf,yf,bf,grdf.x_rho,grdf.y_rho)

h=b_grid * 1.0

hf = h.copy()
smoothing=False
if smoothing:
    print 'smooth the bathymetry and straighten out the edges'
    for i in range(1, h.shape[1]-1):
        for j in range(1, h.shape[0]-1):
            hf[j, i] = 0.5*h[j, i] + 0.125*(h[j+1, i] + 
                                      h[j-1, i] + h[j, i+1] + h[j, i-1])
    
else:
    print 'NO SMOOTHING >>>>'


print 'modifying min/max depth'
#h_min = -1.0
h_min =  0.0
h_max = 50.0

hf[isnan(hf)] = h[isnan(hf)]
hf[isnan(hf)] = h_min
hf[hf<h_min] = h_min
hf[hf>h_max] = h_max



#making buffer area
nll=5
#south
for i in range(nll):
    hf[i, :] = hf[nll,:]

#north
for i in range(1, nll):
    hf[-i, :] = hf[-nll,:]

for i in range(nll):
    hf[:, i] = hf[:,nll]


########
for i in range(1, nll):
    hf[:, -i] = hf[:,-nll]

### case specif
#hf[35:,-1:-2] = h_min
#hf[35:, 0:1] = h_min
   
grdf.h = hf.copy()

if vv:
    fig=pl.figure(62)
    fig.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.pcolor(grdf.lon_rho,grdf.lat_rho,hf)
    pl.plot(grdf.lon_vert,grdf.lat_vert,'-k',alpha=0.1)
    pl.plot(grdf.lon_vert.T,grdf.lat_vert.T,'-k',alpha=0.1)
    pl.clim(1, 10)
    pl.draw()
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('Lon. (Deg.)')
    pl.gca().set_ylabel('Lat. (Deg.)')
    pl.savefig('pic/fig6.png',dpi=450)
    #pl.close()    

if True:
    fig=pl.figure(66)
    fig.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.pcolor(grdf.x_rho,grdf.y_rho,hf)
    #pl.plot(grdf.x_vert,grdf.y_vert,'-k',alpha=0.1)
    #pl.plot(grdf.x_vert.T,grdf.y_vert.T,'-k',alpha=0.1)
    pl.clim(0, 7)
    pl.draw()
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('x (m)')
    pl.gca().set_ylabel('y (m)')
    pl.savefig('pic/fig66.png',dpi=450)
    xb1 = [grdf.x_rho.min()+100, grdf.x_rho.min()+100, 
           grdf.x_rho.max()-100, grdf.x_rho.max()-100]
    yb1 = [grdf.y_rho.max()-100, grdf.y_rho.min()+100, 
           grdf.y_rho.min()+100, grdf.y_rho.max()-100]
    betab1 = [0.0, 0.0, 0.0, 0.0]
    grid_tmp = BoundaryInteractor(xb1,yb1,betab1)
    #pl.close()    

###  Create mask for grid
hf=ma.masked_where(hf<=h_min,hf)
grdf.mask[hf.mask]=0
grdf.mask[~hf.mask]=1

### mask extra points#########################################
if vv:
    pl.figure()
    pl.pcolor(grdf.x_rho,grdf.y_rho,grdf.mask)
    xb1 = [grdf.x_rho.min()+100, grdf.x_rho.min()+100, 
           grdf.x_rho.max()-100, grdf.x_rho.max()-100]
    yb1 = [grdf.y_rho.max()-100, grdf.y_rho.min()+100, 
           grdf.y_rho.min()+100, grdf.y_rho.max()-100]
    betab1 = [0.0, 0.0, 0.0, 0.0]
    grid_tmp = BoundaryInteractor(xb1,yb1,betab1)


#sys.exit
#################################################
#
if (h_min==0):
    verts1=[(-1969.6608235660206, 2483.455882352941),
     (-1913.6500363032883, 243.63687782805437),
     (-606.73166683954014, 211.8212669683262),
     (-522.7154859454422, 269.08936651583736),
     (-690.74784773363831, 339.08371040724023),
     (-322.01016492065219, 1446.2669683257923),
     (284.77336375894492, 1522.6244343891403),
     (345.45171662690473, 2120.7579185520367),
     (1306.9702313038051, 2171.6628959276018),
     (1101.5973446737876, 669.96606334841636),
     (793.53801472876148, 689.05542986425371),
     (760.86505549216781, 778.13914027149349),
     (555.49216886215027, 765.41289592760222),
     (975.57307333264043, 71.832579185520444),
     (1246.2918784358458, 46.38009049773791),
     (2310.4968364277547, 886.31221719457017),
     (2226.4806555336563, 2394.3721719457012)]
    
    grdf.mask_polygon(verts1)

#     #top east bou masking to check for 
    verts1=[(1330.2797379032259, 146.79939516129002),
     (2289.2598575905631, 146.17919921874972),
     (2246.0433467741932, 277.21774193548356),
     (2192.1748991935483, 401.96572580645125)]
    grdf.mask_polygon(verts1)

    #top west bou masking to check for
    verts1=[(-1929.0847381369501, 246.44011339709539),
      (-1974.470505729302, 76.243484925776443),
      (-1427.0046841465594, 79.080095400298433),
      (-1418.4948527229935, 164.17840963595791)]
    grdf.mask_polygon(verts1)



##east bou masking
#verts1=[(-1889.1544117647054, 250.00000000000045),
# (-1946.4225113122166, 110.0113122171947),
# (-1399.1940045248864, 141.82692307692332),
# (-1583.7245475113116, 218.18438914027183)]
#grdf.mask_polygon(verts1)
#
###west bou masking
#verts1=[(1880.4098690257351, 232.13686015271449),
# (1853.0186165511877, 171.58777573529366),
# (2297.0452356122732, 193.21244874151537),
# (2263.8874036694001, 252.31988829185477)]
#grdf.mask_polygon(verts1)


grdf.h=grdf.h *  grdf.mask
grdf.h[grdf.mask==0]=-10.0




if vv==True:
    fig=pl.figure(17)
    fig.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.pcolor(grdf.lon_rho,grdf.lat_rho,hf)
    pl.plot(grdf.lon_vert,grdf.lat_vert,'-k',alpha=0.2)
    pl.plot(grdf.lon_vert.T,grdf.lat_vert.T,'-k',alpha=0.2)
    pl.clim(0, 2)
    pl.draw()
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('Lon. (Deg.)')
    pl.gca().set_ylabel('Lat. (Deg.)')
    pl.savefig('pic/fig17.png',dpi=450)
    #pl.close() 


#############################################

if vv :
    fig=pl.figure(16)
    fig.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    img3=pl.pcolor(grdf.lon_rho,grdf.lat_rho,grdf.mask_rho)
    #gca().xlabel('Lon')
    #gca().ylabel('Lat')
    pl.savefig('pic/mask.png',dpi=450)
    pl.close()    

if False:
    print 'pickling grid..'
    grdf.focus = None  # can't pickle functions...
    fg = open('grdf.pickle', 'wb')
    cPickle.dump(grdf, fg)
    fg.close()
if True:
    print 'Writing netcdf file..'
    rom_grd_name='nri_regional_grd_roms_v2.nc'
    roms.write_grd(grdf, rom_grd_name, verbose=True)
    

if vv:
    fig=pl.figure(676)
    fig.clf()
    ax = pl.subplot(111)
    ax.set_aspect(aspect=1)
    #img3=pl.pcolor(grdf.x_rho,grdf.y_rho,hf)
    img3=pl.contourf(grdf.x_rho,grdf.y_rho,hf)
    pl.plot(grdf.x_vert,grdf.y_vert,'-k',alpha=0.1)
    pl.plot(grdf.x_vert.T,grdf.y_vert.T,'-k',alpha=0.1)
    pl.draw()
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    pl.gca().set_xlabel('x (m)')
    pl.gca().set_ylabel('y (m)')
    pl.savefig('pic/fig676.png',dpi=600)



if False:
    print 'Geo bathy image ... '
    #####################################
    lonmin=x.min()
    lonmax=x.max()
    lonmean= (lonmin + lonmax)/2
    
    latmin=y.min()
    latmax=y.max()
    latmean= (latmin + latmax)/2
    
    #   llcrnrlon        longitude of lower left hand corner of the desired map domain (degrees).
    #   llcrnrlat        latitude of lower left hand corner of the desired map  domain (degrees).
    #   urcrnrlon        longitude of upper right hand corner of the desired map domain (degrees).
    #   urcrnrlat        latitude of upper right hand corner of the desired map  domain (degrees).
    
    #m=Basemap(projection='merc',lon_0=13.5,lat_ts=54.2,resolution='i')
    #boundary resolution must be one of 'c','l','i','h' or 'f'
    dlon=0.05
    dlat=0.05
    fig=pl.figure(11)
    fig.clf()
    m = Basemap(projection='merc',
        resolution='f',
        llcrnrlon=lonmin ,
        llcrnrlat=latmin ,
        urcrnrlon=lonmax ,
        urcrnrlat=latmax ,
        lon_0=lonmean,
        lat_0=latmean)
    #m.drawcoastlines()
    # fill continents, set lake color same as ocean color. 
    
    #m.fillcontinents(color='coral',lake_color='aqua')
    # the continents will be drawn on top.
    
    #m.drawmapboundary(fill_color='aqua')
    
    
    m.drawmeridians(arange(lonmin, lonmax, dlon), labels=[0,0,1,0],
           dashes=(None, None), linewidth=0.1, color='0.6',fmt='%1.1f')
    
    m.drawparallels(arange(latmin,latmax, dlat), labels=[0,1,0,0],
           dashes=(None, None), linewidth=0.1, color='0.6',fmt='%1.1f')
    
    #m.bluemarble(scale=0.8)
    
    lonp,latp=m(grdf.lon_rho,grdf.lat_rho)
    lonv,latv=m(grdf.lon_vert,grdf.lat_vert)
    img3=m.pcolor(lonp,latp,grdf.h)
    cb=pl.colorbar(img3)
    color_lab='Depth (m)'
    cb.set_label(color_lab)
    m.plot(lonv,latv,'-k',alpha=0.25)
    m.plot(lonv.T,latv.T,'-k',alpha=0.25)
    pl.savefig('pic/fig_coast_line.png',dpi=450)
    #pl.show()
    pl.close()
    

if True:
    print('GETM Writing NetCDF file')
    #________NETCD writing________________
    b=grdf.h *1.0
    b[grdf.mask==0]=-10.0
    
    nx,ny = shape(grdf.h)
    
    
    #____the rotation of the grid_________
    convx = zeros(grdf.angle.shape)
    convx = grdf.angle/2/pi*360
    
    
    missing_value=-10
    nc = netCDF4.Dataset('topo.nc', 'w', format='NETCDF3_CLASSIC')
    nc.createDimension('x', nx+1)
    nc.createDimension('y', ny+1)
    
    nc.createDimension('x_T', nx)
    nc.createDimension('y_T', ny)
    
    grid_type = nc.createVariable('grid_type', 'i')
    grid_type.long_name = 'horizontal grid type'
    grid_type.option_1 = 'cartesian'
    grid_type.option_2 = 'spherical'
    grid_type.option_3 = 'curvilinear'
    grid_type.option_4 = 'spherical curvilinear'
    grid_type[:] = 3
    
    xx_nc = nc.createVariable('xx', 'float', ('x','y',))
    xx_nc.long_name = 'X Position of Vertices'
    xx_nc[:] = grdf.x_vert
     
    yx_nc = nc.createVariable('yx', 'float', ('x','y',))
    yx_nc.long_name = 'Y Positions of Vertices'
    yx_nc[:] = grdf.y_vert
        
    latx_nc = nc.createVariable('latx', 'float', ('x','y',))
    latx_nc.long_name = 'Latitude Positions of Vertices'
    latx_nc[:] = grdf.lat_vert
        
    lonx_nc = nc.createVariable('lonx', 'float', ('x','y',))
    lonx_nc.long_name = 'Longitude Positions of Vertices'
    lonx_nc[:] = grdf.lon_vert
    
    bathy_nc = nc.createVariable('bathymetry', 'float', ('x_T','y_T',))
    bathy_nc.long_name = 'bathymetry at T-points'
    bathy_nc.units = 'm'
    bathy_nc.valid_min = -4.0
    bathy_nc.valid_max = b_grid.max()
    bathy_nc.missing_value = float(missing_value)
    bathy_nc[:] = b
    
    convx_nc = nc.createVariable('convx', 'float', ('x','y',))
    convx_nc[:] = convx
    
    nc.Description = proj_str
    nc.x_orig=xu_orig
    nc.y_orig=yu_orig
    nc.Created = datetime.datetime.now().isoformat()
    nc.close()

"""