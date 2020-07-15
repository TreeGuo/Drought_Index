import numpy as np
from timeit import default_timer as timer
from datetime import datetime
from netCDF4 import Dataset
from cftime import utime

#invocate github module (download)

from climate_indices import indices
from climate_indices import compute

import matplotlib.path as mpath
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

import cartopy.crs as ccrs
import cartopy.feature as cfeature

#load precipitation

all_file = r'E:\disorder\bujidao\Desktop\datas\datas\index_parameters_1948-2014.nc'
fh = Dataset(all_file,'r')                                  # file handle, open in read only mode
# fh.set_auto_mask(False)

#load lat&lon
lons = fh.variables['lon'][:]
lats = fh.variables['lat'][:]
#load time
nctime = fh.variables['time'][:]
t_unit = fh.variables['time'].units
# print(nctime)
prcp = fh.variables['pre'][:,:,:]
#  undef =-999.0
undef = fh.variables['pre'].missing_value

#loading pet 

pet = fh.variables['pet'][:,:,:]
#  undef =-999.0
undef = fh.variables['pet'].missing_value

#loading awc

awc = fh.variables['awc'][:,:]

try :
    t_cal = fh.variables['time'].calendar
except AttributeError : # Attribute doesn't exist
    t_cal = u"gregorian" # or standard

fh.close()                                                        # close the file

mask = (prcp==undef) #被掩码的值为数据中的缺测值
prcp = np.ma.MaskedArray(prcp, mask=mask) #
prcp = prcp.astype(np.float64)

pet = np.ma.MaskedArray(pet, mask=mask) #掩码数组
pet = pet.astype(np.float64)    #转变数组类型

nt,nlat,nlon = prcp.shape

data_utime   = utime(t_unit, calendar=t_cal)
datevar = data_utime.num2date(nctime)

#calculateing scpdsi

scpdsi = np.zeros((804,36,80))
scpdsi[:,:,:] = np.nan

final_pdsi = np.zeros((804,36,80))
final_pdsi[:,:,:] = np.nan

phdi = np.zeros((804,36,80))
phdi[:,:,:] = np.nan

pmdi = np.zeros((804,36,80))
pmdi[:,:,:] = np.nan

zindex = np.zeros((804,36,80))
zindex[:,:,:] = np.nan

for ilat in np.arange(nlat):
    lat = lats[ilat]
    for ilon in np.arange(nlon):
        one_pr  = prcp[:,ilat,ilon]*0.0393701
        one_pet = pet[:,ilat,ilon]*0.0393701
        one_awc = awc[ilat,ilon]*0.0393701
        if(not np.ma.is_masked(one_pr)):
            scpdsi[:,ilat,ilon], final_pdsi[:,ilat,ilon], phdi[:,ilat,ilon], pmdi[:,ilat,ilon], zindex[:,ilat,ilon] = indices.scpdsi(precip_time_series=one_pr,
                                                                                                                                       pet_time_series=one_pet,
                                                                                                                                       awc=one_awc,
                                                                                                                                       data_start_year=1948,
                                                                                                                                       calibration_start_year=1948,
                                                                                                                                       calibration_end_year=2014,
                                                                                                                                       )

#output scpdsi in nc file

f_w = Dataset('scpdsi.nc','w')   

f_w.createDimension('time',804)   
f_w.createDimension('lat',36)   
f_w.createDimension('lon',80)

#created varible
f_w.createVariable('time',np.int,('time'))  
f_w.createVariable('lat',np.float64,('lat'))  
f_w.createVariable('lon',np.float64,('lon'))

f_w.createVariable( 'scpdsi', np.float64, ('time','lat','lon'))
f_w.createVariable( 'final_pdsi', np.float64, ('time','lat','lon'))
f_w.createVariable( 'phdi', np.float64, ('time','lat','lon'))
f_w.createVariable( 'pmdi', np.float64, ('time','lat','lon'))
f_w.createVariable( 'zindex', np.float64, ('time','lat','lon'))

f_w.variables['time'][:]           = nctime
f_w.variables['lat'][:]            = lats
f_w.variables['lon'][:]            = lons
f_w.variables['scpdsi'][:,:,:]     = scpdsi
f_w.variables['final_pdsi'][:,:,:] = final_pdsi
f_w.variables['phdi'][:,:,:]       = phdi
f_w.variables['pmdi'][:,:,:]       = pmdi
f_w.variables['zindex'][:,:,:]     = zindex

# time.units = "days since 1948-01-01"
# t1.calendar = "standard"

#关闭文件
f_w.close()

