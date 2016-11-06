"""
Function regrids Piomas SIT onto selected grid. Also, a test plotter
is provided to check the regridding
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 21 October 2016
    
Usage
-----
    sitnew,lats,lons = regridPiomas(directory,sit,lat1,lon1,lat2,lon2)
"""

print '\n>>> Using regridPiomas function!'

import numpy as np
import read_SeaIceThick_PIOMAS as CT
from netCDF4 import Dataset
from scipy.interpolate import griddata as g
import matplotlib.colors as c
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm

### Define directories
directorydata1 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/' 
directorydata2 = '/home/zlabe/Surtsey3/'

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)

### Call functions
#lat1,lon1,sit = CT.readPiomas(directorydata1,years,0.0)

### Read new lat/lon grid
files = 'CESM_large_ensemble/SIT/interp_1deg/'
filename = directorydata2 + files + 'b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc'
data = Dataset(filename)
lat2 = data.variables['lat'][:]
lon2 = data.variables['lon'][:]
data.close()

if lat1.ndim == 1:
    if lon1.ndim == 1:
        lon1,lat1 = np.meshgrid(lon1,lat1)
        print 'Made meshgrid of new lats/lons!'

if lat2.ndim == 1:
    if lon2.ndim == 1:
        lon2,lat2 = np.meshgrid(lon2,lat2)
        print 'Made meshgrid of new lats/lons!'

def regrid(lat1,lon1,lat2,lon2,var,years):
    """
    Interpolated on selected grid. Reads PIOMAS in as 4d with 
    [year,month,lat,lon]
    """
    
    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(120*360)))   
    
    varn = np.empty((var.shape[0],var.shape[1],lat2.shape[0],lon2.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(varn.shape[0]):
        for j in xrange(varn.shape[1]):
            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],(lat2,lon2),method='linear')
            varn[i,j,:,:] = z
        print 'Completed: Year %s Regridding---' % (years[i])
    return varn

def netcdfPiomas(lats,lons,var,directory):
    print '\n>>> Using netcdfPIOMAS function!'
    
    name = 'Thickness/piomas_regrid_sit_LENS_19792015.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'PIOMAS SIT from 1979-2015 ' \
                        'interpolated on selected grid from LENS'
    
    ### Dimensions
    ncfile.createDimension('years',var.shape[0])
    ncfile.createDimension('months',var.shape[1])
    ncfile.createDimension('lat',var.shape[2])
    ncfile.createDimension('lon',var.shape[3])
    
    ### Variables
    years = ncfile.createVariable('years','f4',('years'))
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable('sit','f4',('years','months','lat','lon'))
    
    ### Units
    varns.units = 'meters'
    ncfile.title = 'PIOMAS SIT on LENS Grid'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'PIOMAS, University of Washington'
    ncfile.references = 'Zhang and Rothrock [2003]'
    
    ### Data
    years[:] = list(xrange(var.shape[0]))
    months[:] = list(xrange(var.shape[1]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print '*Completed: Created netCDF4 File!'
    
#sitn = regrid(lat1,lon1,lat2,lon2,sit,years)
#netcdfPiomas(lat2,lon2,sitn,directorydata1)
       
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit

fig = plt.figure()
ax = plt.subplot(111)

cmap = ncm.cmap('MPL_cubehelix')

var = sitn[0,0,:,:]

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
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
barlim = np.arange(0,8,1)
values = np.arange(0,7.1,.25)

cs = m.contourf(lon2,lat2,var[:,:],
                values,extend='max',latlon=True)
cs1 = m.contour(lon2,lat2,var[:,:],
                values,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)                       
        
cs.set_cmap(cmap)

cbar = m.colorbar(cs,location='right',pad='10%',drawedges=True)
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))  
cbar.ax.tick_params(axis='x', size=.1)
cbar.set_label(r'\textbf{thickness (m)}')

fig.suptitle(r'\textbf{Testing PIOMAS Regrid}')

plt.savefig('/home/zlabe/Desktop/' + 'testing_regrid_PIOMAS_new.png',dpi=300)
'Completed: Script done!'