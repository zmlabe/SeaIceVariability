"""
Script regrids LENS variables onto 1x1 grid for comparisons with sea ice 
thickness products
 
Notes
-----
    Author : Zachary Labe
    Date   : 13 March 2017 
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_var_LENS as LV
import read_SeaIceThick_LENS as lens
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
from scipy.interpolate import griddata as g
from netCDF4 import Dataset

### Define directories
directorydataL = '/home/zlabe/Surtsey3/CESM_large_ensemble/'
directorydataSIT = '/home/zlabe/Surtsey3/' 
directoryfigure = '/home/zlabe/Desktop/'
directorydata2 = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot LENS AD - %s----' % titletime 

ense = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 1920
year2 = 2080
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
yearslens = np.arange(year1,year2+1,1)          
yearsclimo = np.arange(1981,2010+1,1)

### Select variable
variable = 'SIC'
          
### Read in functions
var,lats1,lons1 = LV.readLENSEnsemble(directorydataL,variable) 
sith,lats2,lons2 = lens.readLENSEnsemble(directorydataSIT,0.15,'historical')
          
### 2D lat/lon arrays          
lons2,lats2 = np.meshgrid(lons2,lats2)
lons1,lats1 = np.meshgrid(lons1,lats1)
          
##########################################################################  
##########################################################################
##########################################################################
### Regrid
def regrid(lat1,lon1,lat2,lon2,var,years):
    """
    Interpolated on selected grid. 
    [year,month,lat,lon]
    """
    
    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(14*144)))   
    
    varn = np.empty((var.shape[0],var.shape[1],lat2.shape[0],lon2.shape[1]))
    
    print 'Completed: Start regridding process:'
    
    for i in xrange(varn.shape[0]):
        for j in xrange(varn.shape[1]):
            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],(lat2,lon2),method='linear')
            varn[i,j,:,:] = z
        print 'Completed: Year %s Regridding---' % (years[i])
    return varn
    
varn = np.empty((var.shape[0],var.shape[1],var.shape[2],sith.shape[3],
                 sith.shape[4]))
for i in xrange(len(ense)):
    varn[i,:,:,:,:] = regrid(lats1,lons1,lats2,lons2,var[i],yearslens)
    print 'Completed: Regridded Ensemble #%s for %s!' % (ense[i], variable)

def netcdfLENS(lats,lons,var,directory):
    print '\n>>> Using netcdf4LENS function!'
    
    name = 'lens_regrid_SIC_19202080.nc'
    filename = directory + name
    ncfile = Dataset(filename,'w',format='NETCDF4')
    ncfile.description = 'LENS SIC interpolated on 1x1 grid' 
    
    ### Dimensions
    ncfile.createDimension('ensemble',var.shape[0])
    ncfile.createDimension('years',var.shape[1])
    ncfile.createDimension('months',var.shape[2])
    ncfile.createDimension('lat',var.shape[3])
    ncfile.createDimension('lon',var.shape[4])
    
    ### Variables
    ensemble = ncfile.createVariable('ensemble','f4',('ensemble'))
    years = ncfile.createVariable('years','f4',('years'))
    months = ncfile.createVariable('months','f4',('months'))
    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
    varns = ncfile.createVariable('lhflx','f4',('ensemble','years','months','lat','lon'))
    
    ### Units
    varns.units = 'fraction of water covered by sea ice'
    ncfile.title = 'LENS SIC'
    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
    ncfile.source = 'NCAR LENS'
    ncfile.references = 'Kay et al. [2013]'
    
    ### Data
    ensemble[:] = list(xrange(var.shape[0]))
    years[:] = list(xrange(var.shape[1]))
    months[:] = list(xrange(var.shape[2]))
    latitude[:] = lats
    longitude[:] = lons
    varns[:] = var
    
    ncfile.close()
    print '*Completed: Created netCDF4 File!'

netcdfLENS(lats2,lons2,varn,directorydata2)