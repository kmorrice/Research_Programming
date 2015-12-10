import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pyproj
from skimage.filters import roberts, sobel, scharr, prewitt
#from scipy import stats


LL_WGS84 = pyproj.Proj(proj='latlong',datum='WGS84',errcheck=True)
SPCS_N_OR = pyproj.Proj(init = 'nad27:3601')

def convertCoords(x, y):
    """Converts model coordinaates x, y from spcs to lon/lat"""
    lon_mod, lat_mod = pyproj.transform(SPCS_N_OR, LL_WGS84, x, y)
    return lon_mod, lat_mod
    
def getMonthSatRegion(filen, lonmin, lonmax, latmin, latmax):
    '''Extracts a particular region of interest from the satellite data.  
    Arguments include the filename, min and max longitude, and min and max 
    latitude.
    Assumes a monthly dataset, with only one record.
    Longitude inputs should be in range of -155:-105, and latitude inputs
    should be in range of 22:51'''

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
    
    return lon_sat, lat_sat, sst_sat
  
def readSELFEslab(fn):
    fn = nc.Dataset(fn,'r')
    x = fn.variables['x'][:]
    y = fn.variables['y'][:]
    conn = fn.variables['connectivity'][:]
    temp = fn.variables['water_temperature'][:]
    [lon_mod, lat_mod] = convertCoords(x, y)
    triang = tri.Triangulation(lon_mod, lat_mod, conn)
    
    return lon_mod, lat_mod, temp, triang

def month_avg(data):
    data_mean = np.nanmean(data, axis = 1)
    return data_mean
    
def interpolateSELFE(triang,data,lon_sat,lat_sat):
    """Interpolates the SELFE data set to the same grid as the satellite data
        grid
    """
    interp = tri.LinearTriInterpolator(triang, data)
    xx,yy = np.meshgrid(lon_sat, lat_sat)
    mod_interp = interp(xx,yy)
    return mod_interp
    
def load_coast(file_name):
    data = np.loadtxt(file_name)
    coast_x = data[:,0]
    coast_y = data[:,1]
    return coast_x, coast_y

def plume_region(lon, lat, data):
    idx_lon = np.array(np.where((lon > -125.2) & (lon < -123.3)))
    idx_lat = np.array(np.where((lat > 45) & (lat < 47)))
    lon_min = idx_lon.min()
    lon_max = idx_lon.max()
    lat_min = idx_lat.min()
    lat_max = idx_lat.max()
    lon_plume = lon[lon_min:lon_max]
    lat_plume = lat[lat_min:lat_max]
    data_plume = data[lat_min:lat_max, lon_min:lon_max]

    return lon_plume, lat_plume, data_plume

def plume_edges(data,type):
    ''' Where type is a type of filter for the canny edge
    '''    
    data = data.astype(float)    
    edge_data = type(np.flipud(data))
    return edge_data

    # Need to incorporate if/else statement for satellite data
    # if data_plume is masked
    # data_plume[np.where(data.mask ==True)] = np.nan
    # data_plume = data_plume.astype(float)
    # else
    # if data_plume is not masked

    
# Read in SELFE data for July 2005
[lon_mod,lat_mod,temp,triang]=readSELFEslab('SELFE_slabs/slab_temp_s-1_2005-07-01_2005-08-01.nc')

# Read in satellite data for July 2005
sat_file = 'SST_sat_data/MW2005182_2005212_sstd.nc'
[lon_sat, lat_sat, sst_sat] = getMonthSatRegion(sat_file,
min(lon_mod),max(lon_mod),min(lat_mod),max(lat_mod))

# Read in coastline for model domain
[coast_x, coast_y] = load_coast('coastline.dat')

# Calculate monthly average for model data and interpolate to satellite grid
temp_mean = month_avg(temp)
temp_mod_interp = interpolateSELFE(triang, temp_mean, lon_sat, lat_sat)


# Create some plots *****************************

plt.figure()#figsize = (16,8))
plt.subplot(1,2,1)
p1 = plt.pcolor(lon_sat, lat_sat, sst_sat, vmin = 10, vmax = 18, cmap = 'Spectral_r')
plt.colorbar(p1,label='Temperature $^o$C')
plt.title('Satellite, July 2005')
plt.subplot(1,2,2)
p2 = plt.pcolor(lon_sat,lat_sat,temp_mod_interp, vmin = 10, vmax = 18, cmap = 'Spectral_r')
plt.colorbar(p2, label='Temperature $^o$C')
plt.title('SELFE Model, July 2005')
plt.tight_layout()
#plt.savefig('Images/July_sat_mod.png',dpi=300)

# Calculate and plot difference between satellite and model data

sat_mod_dif = sst_sat - temp_mod_interp
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

# Plume behavior for satellite
[lon_plume,lat_plume,data_plume] = plume_region(lon_sat,lat_sat,sst_sat)
data_plume[np.where(data_plume.mask == True)] = np.nan
sat_edge = plume_edges(data_plume, sobel)
plt.imshow(sat_edge)
#edge_sobel = sobel(data) # Make sure it's float
#plt.imshow(edge_sobel)

# Plume behavior for model
[lon_plume,lat_plume,data_plume] = plume_region(lon_sat,lat_sat,temp_mod_interp)
edge_data = plume_edges(data_plume,sobel)
plt.imshow(edge_data)

  # Using tripcolor on original model requires that I mask values
#def mask_slab(data):
#    idx_mask = np.isnan(data)
#    return idx_mask  
# Convert model data to correct lat/lon coordinates
#plt.tripcolor(triang, temp_mean)