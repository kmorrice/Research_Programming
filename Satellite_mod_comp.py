# Satellite_mod_comp.py
# Katie Morrice
# Research Programming Project

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pyproj
from skimage.filters import roberts, sobel, scharr, prewitt
import os
import pandas as pd

plt.rcParams.update({'font.size': 14})
    
LL_WGS84 = pyproj.Proj(proj='latlong',datum='WGS84',errcheck=True)
SPCS_N_OR = pyproj.Proj(init = 'nad27:3601')

def convertCoords(x, y):
    '''Converts model coordinaates x, y from state plane coordinates to lon/lat
    '''
    lon_mod, lat_mod = pyproj.transform(SPCS_N_OR, LL_WGS84, x, y)
    return lon_mod, lat_mod
    
def getMonthSatRegion(filen):
    '''Extracts a particular region of interest from the satellite data.
    Assumes a monthly dataset, with only one record.
    Longitude ranges from -155 to -105, and latitude ranges from 22 to 51 for 
    whole data.
    '''
    fn = nc.Dataset(filen,'r')
    
    # Read in the variables of interest. 
    lon = fn.variables['lon'][:]
    lat = fn.variables['lat'][:]
    MWsstd = fn.variables['MWsstd'][0]
    
    # Make correction to longitude
    lon = lon - 360.0
    idx_lon = np.array(np.where((lon > lonmin) & (lon < lonmax)))
    idx_lat = np.array(np.where((lat > latmin) & (lat < latmax)))
    
    lon_min = idx_lon.min()
    lon_max = idx_lon.max()
    lat_min = idx_lat.min()
    lat_max = idx_lat.max()
    
    lon_sat = lon[lon_min:lon_max]
    lat_sat = lat[lat_min:lat_max]
    sst_sat = MWsstd[0, lat_min:lat_max, lon_min:lon_max]
    
    fn.close()
    return lon_sat, lat_sat, sst_sat

def readSELFEslab(fn):
    '''Reads in variables of interest from model data. 
    Coordinates must be converted from state plane coordinates to latitude and 
    longitude. The triangulation is required because the model is an 
    unstructured triangular grid.
    '''
    fn = nc.Dataset(fn,'r')
    x = fn.variables['x'][:]
    y = fn.variables['y'][:]
    conn = fn.variables['connectivity'][:]
    temp = fn.variables['water_temperature'][:]
    [lon_mod, lat_mod] = convertCoords(x, y)
    triang = tri.Triangulation(lon_mod, lat_mod, conn)
    fn.close()
    return lon_mod, lat_mod, temp, triang
    
def month_avg(data):
    '''Calculates monthly average for the model data.  
    '''
    data_m1 = np.ma.masked_array(data,np.isnan(data))
    data_mean = np.mean(data_m1, axis = 1)
    mask = np.ma.getmask(data_mean)
    data_mean[mask] = np.nan
    return data_mean
    
def interpolateSELFE(triang,data,lon_sat,lat_sat):
    '''Interpolates the SELFE data set to the same grid as the satellite. 
    Interpolation method used is due to the model being an unstructured 
    triangular grid.
    '''
    interp = tri.LinearTriInterpolator(triang, data)
    xx,yy = np.meshgrid(lon_sat, lat_sat)
    mod_interp = interp(xx,yy)
    return mod_interp
    
def load_coast(file_name):
    '''Loads coastline data to ease interpretation of plots.
    '''
    data = np.loadtxt(file_name)
    coast_x = data[:,0]
    coast_y = data[:,1]
    return coast_x, coast_y

def plume_region(lon, lat, sat_data, mod_data):
    '''Extracts data for region that is closer to the plume.
    '''
    idx_lon = np.array(np.where((lon > -125.2) & (lon < -123.3)))
    idx_lat = np.array(np.where((lat > 45) & (lat < 47)))
    lon_min = idx_lon.min()
    lon_max = idx_lon.max()
    lat_min = idx_lat.min()
    lat_max = idx_lat.max()
    lon_plume = lon[lon_min:lon_max]
    lat_plume = lat[lat_min:lat_max]
    sat = sat_data[lat_min:lat_max, lon_min:lon_max]
    mod = mod_data[lat_min: lat_max, lon_min:lon_max]

    return lon_plume, lat_plume, sat, mod

def plume_edges(data,type):
    ''' Detects edges in the supplied data set. Type refers to the filter type 
    used for canny edge. 
    '''   
    data1=data.copy()
    data1[np.where(data.mask == True)] = np.nan
    edge_data = type(data1)
    edge_data = np.ma.masked_invalid(edge_data)
    return edge_data

def plot_sat_mod(sat, mod, month, tmin, tmax):
    '''The satellite longitude and latitude are constant across the datasets, so there is no need to 
    supply the lon/lat for each month. tmin and tmax refer to the minimum and maximum used in the plot.
    '''
    # Mask part of satellite domain where there are not model data
    sat_m = np.ma.masked_where(np.ma.getmask(mod), sat)
    plt.figure(figsize=(16,8))
    plt.subplot(1,2,1)
    p1 = plt.pcolor(lon_sat, lat_sat, sat_m, vmin = tmin, vmax = tmax, cmap = 'Spectral_r')
    plt.colorbar(p1,label='Temperature $^o$C')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Satellite, ' + month + ' 2005')

    plt.subplot(1,2,2)
    p2 = plt.pcolor(lon_sat, lat_sat, mod, vmin = tmin, vmax = tmax, cmap = 'Spectral_r')
    plt.colorbar(p2, label= 'Temperature $^o$C')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('SELFE Model, ' + month + ' 2005')
    plt.tight_layout()  
    plt.savefig('Images/Full/sat_mod_' + month + '.png',dpi = 300)
    plt.close()

def plot_sat_mod_plume(sat, mod, month, tmin, tmax):
    '''Plots the plume region off of the Columbia river. ''' 
    plt.figure(figsize=(16,8))
    plt.subplot(1,2,1)
    p1 = plt.pcolor(lon_p, lat_p, sat, vmin = tmin, vmax = tmax, cmap = 'Spectral_r')
    plt.colorbar(p1,label='Temperature $^o$C')
    plt.plot(coast_x,coast_y,'black')
    plt.xlim(min(lon_p),max(lon_p))
    plt.ylim(min(lat_p),max(lat_p))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Plume, Satellite, ' + month + ' 2005')

    plt.subplot(1,2,2)
    p2 = plt.pcolor(lon_p, lat_p, mod, vmin = tmin, vmax = tmax, cmap = 'Spectral_r')
    plt.colorbar(p2, label= 'Temperature $^o$C')
    plt.plot(coast_x,coast_y,'black')
    plt.xlim(min(lon_p),max(lon_p))
    plt.ylim(min(lat_p),max(lat_p))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Plume, SELFE Model, ' + month + ' 2005')
    plt.tight_layout()
    plt.savefig('Images/Plume/sat_mod_' + month + '.png',dpi = 300)
    plt.close()
    
def plot_edges(sat, mod, month):
    '''Plots the output from the edge detection.'''
    plt.figure(figsize=(16,8))
    plt.subplot(1,2,1)
    plt.pcolor(lon_p,lat_p,sat)
    plt.plot(coast_x,coast_y,'black')
    plt.xlim(min(lon_p),max(lon_p))
    plt.ylim(min(lat_p),max(lat_p))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Satellite, ' + month + ' 2005')
    
    plt.subplot(1,2,2)
    plt.pcolor(lon_p,lat_p,mod)
    plt.plot(coast_x,coast_y,'black')
    plt.xlim(min(lon_p),max(lon_p))
    plt.ylim(min(lat_p),max(lat_p))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('SELFE Model, ' + month + ' 2005')
    plt.tight_layout()
    plt.savefig('Images/Plume/plume_edge_' + month + '.png',dpi = 300)
    plt.close()
    
def check_dir(name):
    if not os.path.isdir(name):
        os.mkdir(name)

check_dir("Images")
check_dir("Images/Full")
check_dir("Images/Plume")

# Read in coastline for model domain
[coast_x, coast_y] = load_coast('coastline.dat')


# Reading in of data by Month and Plotting of Data

# February 2005**********************
# Model Data
[lon_mod,lat_mod,temp1,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-02-01_2005-03-01.nc')

# Minimum and maximum longitude and latitudes for model domain. Satellite data
# extraction relies on these bounds. They are the same for all months.
lonmin = lon_mod.min()
lonmax = lon_mod.max()
latmin = lat_mod.min()
latmax = lat_mod.max()

# Satellite Data
[lon_sat, lat_sat, sst_Feb] = getMonthSatRegion('SST_sat_data/MW2005032_2005059_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_Feb = month_avg(temp1)
mod_Feb = interpolateSELFE(triang, temp_Feb, lon_sat, lat_sat)

plot_sat_mod(sst_Feb, mod_Feb, 'February', 8, 18)

# Extract data for plume region and plot
[lon_p, lat_p, sat_p_Feb, mod_p_Feb] = plume_region(lon_sat, lat_sat, sst_Feb, mod_Feb)
plot_sat_mod_plume(sat_p_Feb, mod_p_Feb, 'February', 10, 18)


# March 2005**********************
# Model Data
[lon_mod,lat_mod,temp2,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-03-01_2005-04-01.nc')

# Satellite Data
[lon_sat, lat_sat, sst_Mar] = getMonthSatRegion('SST_sat_data/MW2005060_2005090_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_Mar = month_avg(temp2)
mod_Mar = interpolateSELFE(triang, temp_Mar, lon_sat, lat_sat)
plot_sat_mod(sst_Mar, mod_Mar, 'March', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_Mar, mod_p_Mar] = plume_region(lon_sat,lat_sat,sst_Mar, mod_Mar)
# Plot plume region
plot_sat_mod_plume(sat_p_Mar, mod_p_Mar, 'March', 8, 14)


# April 2005**********************
# Model Data
[lon_mod,lat_mod,temp3,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-04-01_2005-05-01.nc')

# Satellite Data
[lon_sat, lat_sat, sst_Apr] = getMonthSatRegion('SST_sat_data/MW2005091_2005120_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_Apr = month_avg(temp3)
mod_Apr = interpolateSELFE(triang, temp_Apr, lon_sat, lat_sat)

plot_sat_mod(sst_Apr, mod_Apr, 'April', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_Apr, mod_p_Apr] = plume_region(lon_sat,lat_sat,sst_Apr, mod_Apr)
# Plot plume region
plot_sat_mod_plume(sat_p_Apr, mod_p_Apr, 'April', 8, 14)

# May 2005**********************
# Model Data
[lon_mod,lat_mod,temp4,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-05-01_2005-06-01.nc')

# Satellite Data
[lon_sat, lat_sat, sst_May] = getMonthSatRegion('SST_sat_data/MW2005121_2005151_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_May = month_avg(temp4)
mod_May = interpolateSELFE(triang, temp_May, lon_sat, lat_sat)

plot_sat_mod(sst_May, mod_May, 'May', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_May, mod_p_May] = plume_region(lon_sat,lat_sat,sst_May, mod_May)
# Plot plume region
plot_sat_mod_plume(sat_p_May, mod_p_May, 'May', 10, 16)


# June 2005**********************

# Model Data
[lon_mod,lat_mod,temp5,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-06-01_2005-07-01.nc')

# Satellite Data
[lon_sat, lat_sat, sst_June] = getMonthSatRegion('SST_sat_data/MW2005152_2005181_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_June = month_avg(temp5)
mod_June = interpolateSELFE(triang, temp_June, lon_sat, lat_sat)

plot_sat_mod(sst_June, mod_June, 'June', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_June, mod_p_June] = plume_region(lon_sat,lat_sat,sst_June, mod_June)
# Plot plume region
plot_sat_mod_plume(sat_p_June, mod_p_June, 'June', 10, 18)

# July 2005**********************

# Model Data
[lon_mod,lat_mod,temp6,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-07-01_2005-08-01.nc')

# Satelite Data
[lon_sat, lat_sat, sst_July] = getMonthSatRegion('SST_sat_data/MW2005182_2005212_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_July = month_avg(temp6)
mod_July = interpolateSELFE(triang, temp_July, lon_sat, lat_sat)

plot_sat_mod(sst_July, mod_July, 'July', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_July, mod_p_July] = plume_region(lon_sat,lat_sat,sst_July, mod_July)
# Plot plume region
plot_sat_mod_plume(sat_p_July, mod_p_July, 'July', 10, 18)

# August 2005**********************

# Model Data
[lon_mod,lat_mod,temp7,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-08-01_2005-09-01.nc')

# Satelite Data
[lon_sat, lat_sat, sst_Aug] = getMonthSatRegion('SST_sat_data/MW2005213_2005243_sstd.nc')

# Calculate monthly average for model data and interpolate to satellite grid
temp_Aug = month_avg(temp7)
mod_Aug = interpolateSELFE(triang, temp_Aug, lon_sat, lat_sat)

plot_sat_mod(sst_Aug, mod_Aug, 'August', 8, 18)

# Extract data for plume region
[lon_p,lat_p,sat_p_Aug, mod_p_Aug] = plume_region(lon_sat,lat_sat,sst_Aug, mod_Aug)
# Plot plume region
plot_sat_mod_plume(sat_p_Aug, mod_p_Aug, 'August', 11, 18)


# Calculate and plot difference between satellite and model data

sat_mod_dif = sst_July - mod_July
plt.figure()
plt.pcolor(lon_sat,lat_sat, sat_mod_dif, vmin = -6, vmax = 6, cmap = 'RdBu')
plt.colorbar()

# Plot difference between satellite and model data, with a focus on the estuary
plt.figure()
plt.pcolor(lon_sat,lat_sat, sat_mod_dif, vmin = -3, vmax = 3, cmap = 'RdBu')
plt.colorbar()
plt.xlim(-125,-123.7)
plt.ylim(45.5, 46.8)
plt.plot(coast_x,coast_y,'k',linewidth=2)

# Plot satellite vs. model data for all months

plt.figure(figsize = (16,8))
plt.subplot(2,4,1)
p1 = plt.plot(sst_Feb, mod_Feb,'.b')
plt.ylabel('Model')
plt.xlim(4,15)
plt.ylim(4,15)
plt.plot([4,15],[4,15],'-k')
plt.title('February')

plt.subplot(2,4,2)
p1 = plt.plot(sst_Mar, mod_Mar,'.b')
plt.xlim(4,15)
plt.ylim(4,15)
plt.plot([4,15],[4,15],'-k')
plt.title('March')

plt.subplot(2,4,3)
p1 = plt.plot(sst_Apr, mod_Apr,'.b')
plt.xlim(4,15)
plt.ylim(4,15)
plt.plot([4,15],[4,15],'-k')
plt.title('April')

plt.subplot(2,4,5)
p1 = plt.plot(sst_May, mod_May,'.b')
plt.ylabel('Model')
plt.xlabel('Satellite')
plt.xlim(9,19)
plt.ylim(9,27)
plt.plot([9,27],[9,27],'-k')
plt.title('May')

plt.subplot(2,4,6)
p1 = plt.plot(sst_June, mod_June,'.b')
plt.xlabel('Satellite')
plt.xlim(9,19)
plt.ylim(9,27)
plt.plot([9,27],[9,27],'-k')
plt.title('June')

plt.subplot(2,4,7)
p1 = plt.plot(sst_July, mod_July,'.b')
plt.xlabel('Satellite')
plt.xlim(9,19)
plt.ylim(9,27)
plt.plot([9,27],[9,27],'-k')
plt.title('July')


plt.subplot(2,4,8)
p1 = plt.plot(sst_Aug, mod_Aug,'.b')
plt.xlabel('Satellite')
plt.xlim(9,19)
plt.ylim(9,27)
plt.plot([9,27],[9,27],'-k')
plt.title('August')
plt.savefig('Images/mod_sat.png',dpi=300)

# Edge Detection

# For February
Feb_edge_sat = plume_edges(sat_p_Feb,sobel)
Feb_edge_mod = plume_edges(mod_p_Feb,sobel)

plot_edges(Feb_edge_sat, Feb_edge_mod, 'February')

# For July
July_edge_sat = plume_edges(sat_p_July,sobel)
July_edge_mod = plume_edges(mod_p_July,sobel)

plot_edges(July_edge_sat, July_edge_mod, 'July')

# Statistics
sst_stats = np.zeros((7,2))
mod_stats = np.zeros((7,2))

# satellite data to calculate means for
sst_list=[sst_Feb,sst_Mar,sst_Apr,sst_May,sst_June,sst_July,sst_Aug]

# model data to calculate means for
mod_list=[mod_Feb,mod_Mar,mod_Apr,mod_May,mod_June,mod_July,mod_Aug]
          
def compute_stats(data):
    m = data.mean()
    v = data.var()
    return m, v
    
for i in range(0,7):
    [sst_stats[i,0],sst_stats[i,1]] = compute_stats(sst_list[i])
    [mod_stats[i,0],mod_stats[i,1]] = compute_stats(mod_list[i])         
          
# Creation and concatenation of satellite and model statistics
df1 = pd.DataFrame(sst_stats, index=['February','March','April','May','June','July','August'])
df1.columns = ['Sat. mean','Sat. variance']
df2 = pd.DataFrame(mod_stats, index = ['February','March','April','May','June','July','August'])
df2.columns = ['Mod 1 mean','Mod 1 variance']
stats_df = pd.concat([df1, df2], axis = 1)

# Read in additional model run statistics
stats_run2 = pd.read_pickle('run19.pkl')
stats_run2.drop(stats_run2.columns[[0,1]], axis = 1, inplace = True)
stats_run2

# Concatenate the two data frames
stats_df2 = pd.concat([stats_df, stats_run2], axis = 1)

# Reorder columns
stats_df3 = stats_df2.reindex(columns=['Sat. mean','Mod 1 mean','Mod 2 mean','Sat. variance',
                                       'Mod 1 variance','Mod 2 variance'])
print(stats_df3)

dates = np.arange('2005-02', '2005-09', dtype = 'datetime64[M]')
sat = stats_df3['Sat. mean'].values
run16 = stats_df3['Mod 1 mean'].values
run19 = stats_df3['Mod 2 mean'].values

plt.figure(figsize = (8,8))
plt.plot(dates,sat,'or', label="Satellite")
plt.plot(dates,run16,'ob', label = "run16")
plt.plot(dates,run19,'og', label = "run19")
plt.legend(loc = 'upper left')
plt.ylabel('Mean Temperature $^o$C')
plt.xlim(731970,732167)
plt.xticks(rotation = 45)
plt.tight_layout()
plt.savefig('Images/monthly_means.png',dpi=300)