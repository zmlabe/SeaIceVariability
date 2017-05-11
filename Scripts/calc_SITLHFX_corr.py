"""
Scripts calculates spatial correlations between 2 m temperature
and Arctic sea ice thickness in LENS

Notes
-----
    Author : Zachary Labe
    Date   : 16 February 2017
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from scipy import signal
import read_var_LENS as LV
import read_SeaIceThick_LENS as lens
from mpl_toolkits.basemap import Basemap

### Define directories
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
print '\n' '----Plot LENS Correlations - %s----' % titletime 

### Alott time series
year1 = 1920
year2 = 2080
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ense = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))
          
#### Read in functions
#sith,lats,lons = lens.readLENSEnsemble(directorydataSIT,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydataSIT,0.15,'rcp85')

### Read T2M
data = Dataset(directorydata2 + 'lens_regrid_LHFLX_19202080.nc')
lhfxall = data.variables['lhflx'][:]
data.close()

### Combine SIT periods
sitall = np.append(sith,sitf,axis=1)
          
#### 2D lat/lon arrays          
lons,lats = np.meshgrid(lons,lats)

### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

###########################################################################
###########################################################################
###########################################################################
### Call functions
def deTrend(y):
    x = np.arange(y.shape[0])
    
    slopes = np.empty((y.shape[1],y.shape[2]))
    intercepts = np.empty((y.shape[1],y.shape[2]))
    for i in xrange(y.shape[1]):
        for j in xrange(y.shape[2]):
            mask = np.isfinite(y[:,i,j])
            yy = y[:,i,j]           
            
            if np.isfinite(np.nanmean(yy)):
                slopes[i,j], intercepts[i,j], r_value, p_value, std_err = sts.linregress(x[mask],yy[mask])
            else:
                slopes[i,j] = np.nan
                intercepts[i,j] = np.nan
    
    y_detrend = np.empty(y.shape)        
    for i in xrange(y.shape[0]):
        y_detrend[i,:,:] = y[i,:,:] - (slopes*x[i] + intercept)

#    y_detrend = np.empty(y.shape)
#    for i in xrange(y.shape[1]):
#        for j in xrange(y.shape[2]):
#            y_detrend[:,i,j] = signal.detrend(y[:,i,j],type='linear')
     
    print 'Completed: Detrended SIV data!' 
    return y_detrend       
    
sit_w = np.nanmean(sitall[:,:,0:3,:,:],axis=2)
sit_sp = np.nanmean(sitall[:,:,3:6,:,:],axis=2)
sit_su = np.nanmean(sitall[:,:,6:9,:,:],axis=2)
sit_f = np.nanmean(sitall[:,:,9:12,:,:],axis=2)

lhfx_w = np.nanmean(lhfxall[:,:,0:3,:,:],axis=2)
lhfx_sp = np.nanmean(lhfxall[:,:,3:6,:,:],axis=2)
lhfx_su = np.nanmean(lhfxall[:,:,6:9,:,:],axis=2)
lhfx_f = np.nanmean(lhfxall[:,:,9:12,:,:],axis=2)

def corr(sit,tas):
    varx = sit
    vary = tas
    corr = np.empty((sit.shape[1],sit.shape[2]))
    for i in xrange(sit.shape[1]):
        for j in xrange(sit.shape[2]):
            corr[i,j] = sts.stats.pearsonr(varx[:,i,j],vary[:,i,j])[0]
        
    corr[np.where(corr == 1.)] = np.nan
    
    print 'Completed: Correlated SIV and AD data!'
    return corr

### Compute correlations
corrn_w = np.empty((sitall.shape[0],sitall.shape[3],sitall.shape[4]))
corrn_sp = np.empty((sitall.shape[0],sitall.shape[3],sitall.shape[4]))
corrn_su = np.empty((sitall.shape[0],sitall.shape[3],sitall.shape[4]))
corrn_f = np.empty((sitall.shape[0],sitall.shape[3],sitall.shape[4]))
for i in xrange(sitall.shape[0]):
    corrn_w[i,:,:] = corr(sit_w[i],lhfx_w[i])
    corrn_sp[i,:,:] = corr(sit_sp[i],lhfx_sp[i])
    corrn_su[i,:,:] = corr(sit_su[i],lhfx_su[i])
    corrn_f[i,:,:] = corr(sit_f[i],lhfx_f[i])

    print 'Completed: Correlations for ensemble %s!' % ense[i]
    
corr_w = np.nanmean(corrn_w,axis=0)
corr_sp = np.nanmean(corrn_sp,axis=0)
corr_su = np.nanmean(corrn_su,axis=0)
corr_f = np.nanmean(corrn_f,axis=0)

#### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-1,1.1,.5)
values = np.arange(-1,1.1,0.1)

cs = m.contourf(lons,lats,corr_w,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_w,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cs.set_cmap('RdBu_r')
ax.annotate(r'\textbf{JFM}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_sp,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_sp,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cs.set_cmap('RdBu_r')

ax.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(223)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_su,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_su,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cs.set_cmap('RdBu_r')
ax.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(224)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_f,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_f,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cs.set_cmap('RdBu_r')

ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='Both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{Correlation Coefficient}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'sitlhfx_corrs.png',dpi=300)