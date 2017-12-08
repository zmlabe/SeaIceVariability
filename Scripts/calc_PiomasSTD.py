"""
Script calculates and plots STD in PIOMAS SIT

Author : Zachary Labe
Date : 13 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import read_SeaIceThick_PIOMAS as CT
import nclcmaps as ncm
import scipy.stats as sts

### Define directories
directorydata = '/surtsey/zlabe/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate PIOMAS STD - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',r'Sep',
          r'Oct',r'Nov',r'Dec']

### Call functions
lats,lons,sitq = CT.readPiomas(directorydata,years,0.15)

def deTrend(y):
    x = np.arange(y.shape[0])
    
    slopes = np.empty((y.shape[1],y.shape[2],y.shape[3]))
    intercepts = np.empty((y.shape[1],y.shape[2],y.shape[3]))
    for mo in xrange(y.shape[1]):
        for i in xrange(y.shape[2]):
            for j in xrange(y.shape[3]):
                mask = np.isfinite(y[:,mo,i,j])
                yy = y[:,mo,i,j]           
                
                if np.isfinite(np.nanmean(yy)):
                    slopes[mo,i,j], intercepts[mo,i,j], r_value, p_value, std_err = sts.linregress(x[mask],yy[mask])
                else:
                    slopes[mo,i,j] = np.nan
                    intercepts[mo,i,j] = np.nan
        print 'Regressed over month %s!' % (mo)
    
    y_detrend = np.empty(y.shape)        
    for i in xrange(y.shape[0]):
        y_detrend[i,:,:,:] = y[i,:,:,:] - (slopes*x[i] + intercepts)
        print 'Detrended over year %s!' % (i)
     
    print 'Completed: Detrended SIV data!' 
    return y_detrend
    
sitq_dt = deTrend(sitq)

#### Take monthly mean
def monSTD(sitq,months):
    sitstd = np.zeros((sitq.shape[1],sitq.shape[2],sitq.shape[3]))
    for mo in xrange(sitq.shape[1]):
        sit = sitq[:,mo,:,:]
        for i in xrange(0,sit.shape[1]):
            for j in xrange(0,sit.shape[2]):        
                sitstd[mo,i,j] = np.nanstd(sit[:,i,j])
        print 'Completed: Month std %s done!' % (months[mo])
    print 'Completed: Calculated STD!'
    
    return sitstd

sitstd = monSTD(sitq_dt,months)  

sittrend_w = np.nanmean(sitstd[0:3],axis=0)
sittrend_sp = np.nanmean(sitstd[3:6],axis=0)
sittrend_su = np.nanmean(sitstd[6:9],axis=0)
sittrend_f = np.nanmean(sitstd[9:12],axis=0)

###########################################################################
###########################################################################
###########################################################################
### Plot Composites
### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
            
var = sittrend_w 

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(0,1.1,0.5)
values = np.arange(0,1.1,0.1)

cs = m.contourf(lons,lats,var,
                values,latlon=True,extend='max')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = plt.cm.get_cmap('cubehelix_r')        
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JFM}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)

var = sittrend_sp           
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=0)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,var,
                values,latlon=True,extend='max')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = plt.cm.get_cmap('cubehelix_r')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(223)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
            
var = sittrend_su            
                    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,var,
                values,latlon=True,extend='max')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = plt.cm.get_cmap('cubehelix_r')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(224)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
            
var = sittrend_f          
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,var,
                values,latlon=True,extend='max')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = plt.cm.get_cmap('cubehelix_r')        
cs.set_cmap(cmap)
ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='max',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{std. dev. (meters)}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'sit_std_7915.png',dpi=300)

print 'Completed: Script done!'