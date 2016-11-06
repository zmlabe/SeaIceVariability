"""
Script looks at NCEP/NCAR reanalysis trends
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 27 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_NCEP as NP
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey/NCEP/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot NCEP - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in functions
lats,lons,h5 = NP.readNCEP(directorydata,years,'heights','500')     

### calculate climo
def climo(var,years,yearmin,yearmax):
    """
    Calculates climatology based on given years
    """
    yr = np.where((years >= yearmin) & (years <= yearmax))[0]
    
    meanvar = np.nanmean(var[yr,:,:,:],axis=0)
    
    print 'Completed: Calculated mean climatology!'
    return meanvar

### Calculate anomalies    
meanh5 = climo(h5,years,1981,2010)

year1 = 2010
year2 = 2015
yrq = np.where((years >= year1) & (years <= year2))[0]  

#months = 'AMJ'
months = 'JAS'

#anomh5 = h5[yrq,3:6,:,:] - meanh5[3:6,:,:]
anomh5 = h5[yrq,6:9,:,:] - meanh5[6:9,:,:]

summeranom = np.nanmean(np.nanmean(anomh5,axis=1),axis=0)

var = summeranom

### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

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
barlim = np.arange(-100,101,50)

var, lons_cyclic = addcyclic(var, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var[:,:],
                np.arange(-100,101,5),extend='both')
cs1 = m.contour(x,y,var[:,:],
                np.arange(-100,101,5),linewidths=0.2,colors='k',
                linestyles='-')

cmap = ncm.cmap('BlueDarkRed18')         
cs.set_cmap(cmap)

cbar = m.colorbar(cs,location='right',pad='10%',drawedges=True)
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))  
cbar.ax.tick_params(axis='x', size=.1)
cbar.set_label(r'\textbf{geopotential meters}')

fig.suptitle(r'\textbf{500 mb anomalies - %s %s-%s}' % (months,year1,year2))

plt.savefig(directoryfigure + 'H5/h5_%s_%s%s.png' % (months,year1,year2),
            dpi=300)
'Completed: Script done!'