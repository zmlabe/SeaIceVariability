"""
Script calculates ip time series 
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 22 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import calc_PiomasArea as CA
import read_SeaIceProduction as CP
from mpl_toolkits.basemap import Basemap

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate Ice Production - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2013
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

### Call functions
lats,lons,ip = CP.readPiomas(directorydata,years)
area = CA.readPiomasArea(directorydata)



val = ip[-2,1,:,:]




### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

#### Define figure
fig = plt.figure()   
ax = plt.subplot(111)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.1)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[True,True,True,True],
                linewidth=0.2,color='k',fontsize=5)
m.drawmeridians(meridians,labels=[True,True,True,True],
                linewidth=0.2,color='k',fontsize=5)
m.drawlsmask(land_color='darkgrey',ocean_color='w')

### Adjust maximum limits
values = np.arange(-0.2,0.21,0.01)  

### Plot filled contours    
cs = m.contourf(lons[:,:],lats[:,:],val[:,:],
                values,latlon=True,extend='both')
cs1 = m.contour(lons[:,:],lats[:,:],val[:,:],
                values,latlon=True,linestyles='-',
                colors='k',linewidths=0.05)
                  
### Set colormap          
cmap = plt.cm.get_cmap('seismic_r')                      
cs.set_cmap(cmap)
                     
cbar = m.colorbar(cs,drawedges=True,location='bottom',pad = 0.3,
                  extend='both')

cbar.set_label(r'Ice Production [m/day]')
cbar.set_ticks(np.arange(-0.2,0.21,0.1))
cbar.set_ticklabels(map(str,np.arange(-0.2,0.21,0.1))) 

fig.suptitle(r'\textbf{PIOMAS, Ice Production}',fontsize=14)
fig.subplots_adjust(bottom=0.15)

### Save figure
plt.savefig(directoryfigure +'iceproduct.png',dpi=500)

###########################################################################
###########################################################################
###########################################################################

def weightIP(var,area):
    """
    Area weights ip array 4d [year,month,lat,lon] into [year,month]
    """
    ipyr = np.empty((var.shape[0],var.shape[1]))
    for i in xrange(var.shape[0]):
        for j in xrange(var.shape[1]):
            varq = var[i,j,:,:]
            mask = np.isfinite(varq) & np.isfinite(area)
            varmask = varq[mask]
            areamask = area[mask]
            ipyr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted ip average!' 
    return ipyr
    
ipyr = weightIP(ip,area)

plt.figure()
plt.plot(ipyr)