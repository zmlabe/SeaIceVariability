"""
Script plots time series for month of March using satellite data and
modeled data from PIOMAS for sea ice thickness
 
Source 1 : ftp://sidads.colorado.edu/pub/projects/SIPN/seaice_thickness/
Source 2 : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
Author : Zachary Labe
Date : 20 July 2016
"""

### Import Modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import griddata as g
import datetime
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import scipy.stats as sts

### Define directories   
directoryfigure = directoryfigure = '/home/zlabe/Desktop/'
directorydata = '/home/zlabe/Surtsey/seaice_obs/Thk/March/'

### Define constants
years = np.arange(1979,2016,1)
timex = np.arange(0,37,1)

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr

print '\n' 'March SIT Time Series Plot - %s' '\n' % titletime 

def readPIOMAS(directory):
    """
    Read PIOMAS March time series (1979-2015)
    """

    filename = 'piomas_regrid_March_19792015.nc'    
    
    data = Dataset(directory + filename,'r')
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    sitp = data.variables['thick'][:]
    data.close()
    
    print 'Completed: PIOMAS data read!'
    return lats,lons,sitp
    
def readSatG(directory):
    """
    Read satellite data with ICESat-G and CryoSat (2003-2009)
    """
    
    filename = 'satelliteG_regrid_March_20032015.nc'
    
    data = Dataset(directory + filename,'r')
    sitsg = data.variables['thick'][:]
    data.close()
    
    print 'Completed: Satellite-G data read!'
    return sitsg
    
def readSatJ(directory):
    """
    Read satellite data with ICESat-J and CryoSat (2003-2009)
    """
    
    filename = 'satelliteJ_regrid_March_20042015.nc'
    
    data = Dataset(directory + filename,'r')
    sitsj = data.variables['thick'][:]
    data.close()
    
    print 'Completed: Satellite-J data read! \n'
    return sitsj 
    
### Call functions
lats,lons,sitp = readPIOMAS(directorydata)
sitsg = readSatG(directorydata)
sitsj = readSatJ(directorydata)

### Create mask
maskgrid = sitsj.copy()
maskgrid = maskgrid[0,:,:]
maskgrid[np.where(maskgrid > 0.)] = 1.
maskgrid[np.where(maskgrid != 1.)] = 0.

### Complete time series
years19792003 = np.empty((24,180,180))
years19792003.fill(np.nan)
years19792004 = np.empty((25,180,180))
years19792004.fill(np.nan)

sitsg = np.append(years19792003,sitsg,axis=0)
sitsj = np.append(years19792004,sitsj,axis=0)

### Apply mask
sitsg = sitsg * maskgrid
sitsj = sitsj * maskgrid
sitp = sitp * maskgrid

sitsg[np.where(sitsg == 0.)] = np.nan
sitsj[np.where(sitsj == 0.)] = np.nan
sitp[np.where(sitp == 0.)] = np.nan

### Take areal average
mean_sitp = np.nanmean(np.nanmean(sitp,axis=1),axis=1)
mean_sitsg = np.nanmean(np.nanmean(sitsg,axis=1),axis=1)
mean_sitsj = np.nanmean(np.nanmean(sitsj,axis=1),axis=1)

cryo = mean_sitsg[-5:]
years19792010 = np.empty((37-5))
years19792010.fill(np.nan)
mean_cryo = np.append(years19792010,cryo)

### Calculate trends
slope,intercept,r_value,p_value,std_error = sts.linregress(timex,mean_sitp)
linep = slope*timex + intercept

masksg = ~np.isnan(mean_sitsg)
masksg[-6:] = False
slopesg,interceptsg,r_valuesg,p_valuesg,std_errorsg = sts.linregress(timex[masksg],mean_sitsg[masksg])
linesg = slopesg*timex + interceptsg

masksj = ~np.isnan(mean_sitsj)
masksj[-6:] = False
slopesj,interceptsj,r_valuesj,p_valuesj,std_errorsj = sts.linregress(timex[masksj],mean_sitsj[masksj])
linesj = slopesj*timex + interceptsj

maskc = ~np.isnan(mean_cryo)
slopec,interceptc,r_valuec,p_valuec,std_errorc = sts.linregress(timex[maskc],mean_cryo[maskc])
linec = slopec*timex + interceptc

print 'Loss of %s meters per decade [PIOMAS]' % round((slope*10.),2) 
print 'Loss of %s meters per decade [ICESat-G]' % round((slopesg*10.),2) 
print 'Loss of %s meters per decade [ICESat-J]' % round((slopesj*10.),2) 
print 'Loss of %s meters per decade [CryoSat] \n' % round((slopec*10.),2)

### Plot March time series
fig = plt.figure()
ax = plt.subplot(111)
    
### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

### Plot
trendp = plt.plot(linep,linewidth=1.4,linestyle='-',color='seagreen'),
sg = plt.plot(mean_sitsg,linestyle='-',linewidth=0.8,color='saddlebrown',
             label=r'\textbf{ICESat-G}',marker='o',markersize=3)
sj = plt.plot(mean_sitsj,linestyle='-',linewidth=0.8,color='darkslateblue',
             label=r'\textbf{ICESat-J}',marker='o',markersize=3)
c = plt.plot(mean_cryo,linestyle='-',linewidth=0.8,color='fuchsia',
             label=r'\textbf{CryoSat-2}',marker='o',markersize=3)
p = plt.plot(mean_sitp,linestyle='-',linewidth=0.8,color='seagreen',
             label=r'\textbf{PIOMAS}',marker='^',markersize=4)            

### Labels for x/y
labelsy = map(str,np.arange(1,5,1))
labelsx = map(str,np.arange(1979,2016,3))
plt.xticks(np.arange(0,37,3),labelsx)
plt.yticks(np.arange(1,5,1),labelsy)
plt.ylabel(r'\textbf{Thickness (meters)}',fontsize=11)

### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='k',zorder=1,alpha=0.2)

### Add limits to axes
plt.ylim([1,4])
plt.xlim([0,36])

### Add legend
plt.legend(shadow=False,fontsize=9,loc='center',
                       fancybox=True,ncol=4,bbox_to_anchor=(0.5,-0.15),
                        frameon=False)

### Add title
fig.suptitle(r'\textbf{March Average Sea Ice Thickness (1979-2015)}',
             fontsize=14)
             
### Create subplot  
diffg = mean_sitsg - mean_sitp 
diffj = mean_sitsj - mean_sitp  
diffc = mean_cryo[:] - mean_sitp[:] 

zero = [0]*len(diffg)

yearsub = np.arange(2003,2016,2)      
             
a = plt.axes([.18, .18, .25, .25], axisbg='w')
for axis in ['top','bottom','left','right']:
  a.spines[axis].set_linewidth(2)
a.set_axis_bgcolor('darkgrey')
plt.plot(zero,color='k',linewidth=2)        
plt.plot(diffg,color='saddlebrown',marker='o',markersize=3)
plt.plot(diffj,color='darkslateblue',marker='o',markersize=3)
plt.plot(diffc,color='fuchsia',marker='o',markersize=3)

plt.title(r'\textbf{Difference, [Satellite -- PIOMAS]}',fontsize=7)
plt.ylim([-1.,1.])
plt.xlim([24,36])
plt.grid(alpha=0.4)
labelsx2 = map(str,yearsub)
labelsy2 = map(str,np.arange(-1,1.5,0.5))
plt.xticks(np.arange(24,37,2),labelsx2,fontsize=6)
plt.yticks(np.arange(-1,1.5,0.5),labelsy2,fontsize=6)
plt.ylabel(r'\textbf{meters}',labelpad=0.2,fontsize=7)

### Create 2nd subplot
masking = maskgrid.copy()
masking[np.where(masking == 0.)] = np.nan

a2 = plt.axes([.63, .59, .29, .29], axisbg='w')   
m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,resolution='l',round=True)
m.drawmapboundary(fill_color = 'white')
m.drawcoastlines(color = 'darkgrey',linewidth=0.2)
m.drawlsmask(land_color='darkgrey',ocean_color='snow')
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0.25)
m.drawmeridians(meridians,labels=[True,True,True,True],linewidth=0.25,
                fontsize=4)

cs = m.contourf(lons,lats,masking,np.arange(0,3,1),
                latlon=True,colors='seagreen')

fig.subplots_adjust(bottom=0.15)

### Save figure
plt.savefig(directoryfigure + 'March_SIT_timeseries.png',dpi=300)

print 'Completed: Script done!'