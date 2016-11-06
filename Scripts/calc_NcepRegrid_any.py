"""
Function regrids NCEP data onto selected grid. Also, a test plotter
is provided to check the regridding
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 27 October 2016
    
Usage
-----
    varnew,lats,lons = regridNcep(directory,var,lat1,lon1,lat2,lon2)
"""

print '\n>>> Using regridNcep function!'

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata as g
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import read_NCEP as NP
import nclcmaps as ncm

### Define directories
directorydata1 = '/home/zlabe/Surtsey/NCEP/' 
directorydata2 = '/home/zlabe/Surtsey3/'

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)

varnames = 'slp'

### Call functions
lat1,lon1,var = NP.readNCEP(directorydata1,years,varnames,'surface')  

### Read new lat/lon grid
files = 'CESM_large_ensemble/SIT/interp_1deg/'
filename = directorydata2 + files + 'b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc'
data = Dataset(filename)
lat2 = data.variables['lat'][:]
lon2 = data.variables['lon'][:]
data.close()

if lat2.ndim == 1:
    if lon2.ndim == 1:
        lon2,lat2 = np.meshgrid(lon2,lat2)
        print 'Made meshgrid of new lats/lons!'
        
if lat1.ndim == 1:
    if lon1.ndim == 1:
        lon1,lat1 = np.meshgrid(lon1,lat1)
        print 'Made meshgrid of new lats/lons!'

def regrid(lat1,lon1,lat2,lon2,var,years):
    """
    Interpolated on selected grid. Reads NCEP in as 4d with 
    [year,month,lat,lon]
    """
    
    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(73*144)))   
    
    varn = np.empty((var.shape[0],var.shape[1],lat2.shape[0],lon2.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(1):
        for j in xrange(varn.shape[1]):
            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],(lat2,lon2),method='linear')
            varn[i,j,:,:] = np.flipud(np.fliplr(z))
        print 'Completed: Year %s Regridding---' % (years[i])
    return varn

def netcdfNcep(lats,lons,var,varnames):
    print '\n>>> Using netcdfNcep function!'
    
    directory = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'
    name = 'ncep_regrid_%s_LENS_19792015.nc' % varnames
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'NCEP from 1979-2015 ' \
                        'interpolated on selected grid from LENS'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('months',var.shape[1])
    ncfile.createDimension('lat',var.shape[2])
    ncfile.createDimension('lon',var.shape[3])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lat'))
    longitude = ncfile.createVariable('lon','f4',('lon','lon'))
    varns = ncfile.createVariable('sit','f4',('years','months','lat','lon'))
    
    ### Units
    ncfile.title = 'NCEP %s on LENS Grid' % varnames
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NCEP/NCAR Reanalysis 1'
    ncfile.references = 'Kalnay et al. [1996]'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    months[:] = list(xrange(var.shape[1]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print '*Completed: Created netCDF4 File!'
    
varn = regrid(lat1,lon1,lat2,lon2,var,years)
#netcdfNcep(lat2,lon2,varn,directorydata1)
       
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

cmap = ncm.cmap('testcmap') 

fig = plt.figure()
ax = plt.subplot(111)

var = var[0,0,:,:]

m = Basemap(projection='npstere',boundinglat=62,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[True,True,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
cs = m.contourf(lon1,lat1,var[:,:],
                extend='both',latlon=True)
cs1 = m.contour(lon1,lat1,var[:,:],
                linewidths=0.2,colors='k',
                linestyles='-',latlon=True)                       
        
cs.set_cmap(cmap)

cbar = m.colorbar(cs,location='right',pad='10%',drawedges=True) 
cbar.ax.tick_params(axis='x', size=.1)
cbar.set_label(r'\textbf{SLP (mb)}')

fig.suptitle(r'\textbf{Testing NCEP Regrid}')

plt.savefig('/home/zlabe/Desktop/' + 'testing_regrid_ncep.png',dpi=300)
#'Completed: Script done!'