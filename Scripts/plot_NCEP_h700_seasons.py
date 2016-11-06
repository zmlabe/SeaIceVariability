"""
Script looks at NCEP/NCAR reanalysis trends
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 2 November 2016
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
directoryfigure = '/home/zlabe/Desktop/H7/'

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
lats,lons,h7 = NP.readNCEP(directorydata,years,'heights','700')      

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
meanh7 = climo(h7,years,1981,2010)
anomalies = h7 - meanh7

springanom = np.nanmean(anomalies[:,3:6,:,:],axis=1)
summeranom = np.nanmean(anomalies[:,6:9,:,:],axis=1)

#totalanom = np.append(springanom[-6:,:,:],summeranom[-6:,:,:],axis=0)
totalanom = np.append(springanom[-12:-6,:,:],summeranom[-12:-6,:,:],axis=0)

### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for i in xrange(totalanom.shape[0]):
    
    var = totalanom[i,:,:]    
    
    ax = plt.subplot(2,6,i+1)
    
    m = Basemap(projection='npstere',boundinglat=62,lon_0=270,
                resolution='l',round =True)
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.3)
#    parallels = np.arange(50,90,10)
#    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0,color='k',fontsize=6)
#    m.drawmeridians(meridians,labels=[True,True,False,False],
#                    linewidth=0,color='k',fontsize=6)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    # Make the plot continuous
    barlim = np.arange(-100,101,50)
    values = np.arange(-100,101,5)
    
    var, lons_cyclic = addcyclic(var, lons)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
    x, y = m(lon2d, lat2d)
    
    cs = m.contourf(x,y,var[:,:],
                    values,extend='both')
#    cs1 = m.contour(x,y,var[:,:],
#                    values,linewidths=0.2,colors='k',
#                    linestyles='-')
    
    cmap = ncm.cmap('BlueDarkRed18')        
    cs.set_cmap(cmap)

plt.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.01, 0.635),
            xycoords='figure fraction',fontsize=22,color='k')
plt.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(0.03, 0.34),
            xycoords='figure fraction',fontsize=22,color='k')
            
plt.annotate(r'2004', xy=(0, 0), xytext=(0.155, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')
plt.annotate(r'2005', xy=(0, 0), xytext=(0.285, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')
plt.annotate(r'2006', xy=(0, 0), xytext=(0.415, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')
plt.annotate(r'2007', xy=(0, 0), xytext=(0.545, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')
plt.annotate(r'2008', xy=(0, 0), xytext=(0.675, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')
plt.annotate(r'2009', xy=(0, 0), xytext=(0.805, 0.51),
            xycoords='figure fraction',fontsize=13,color='k')            
    
cbar_ax = fig.add_axes([0.312,0.2,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)
cbar.set_label(r'\textbf{H7( m )}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.1)

fig.subplots_adjust(hspace=-0.4)
fig.subplots_adjust(wspace=0)
fig.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure + 'h7_seasons_20042009.png',
            dpi=300)
print 'Completed: Script done!'