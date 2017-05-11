"""
Scripts calculates Arctic Dipole using 2nd eof of MSLP anomaly north
of 70N. 
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 6 October 2016
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from scipy import signal
from mpl_toolkits.basemap import Basemap
import read_SeaIceThick_LENS as lens

### Define directories
directorydatal = '/home/zlabe/Surtsey3/'
directoryfigure = '/home/zlabe/Desktop/'
directorydatap = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----LENS/PIOMAS SIT Correlations - %s----' % titletime 

### Alott time series
yearmin = 1920
yearmax = 2080
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,39,1)) + map(str,np.arange(101,106,1))
          
def readPIOMAS(directorydata,threshold):
    files = 'piomas_regrid_sit_LENS_19792015.nc'
    filename = directorydata + files
    
    data = Dataset(filename)
    sitp = data.variables['sit'][:,:,156:180,:] # lats > 65
    data.close()
    
    ### Mask out threshold values
    if threshold == 'None':
        sitp[np.where(sitp < 0)] = np.nan
        sitp[np.where(sitp > 12)] = np.nan
    else:
        sitp[np.where(sitp < threshold)] = np.nan
        sitp[np.where(sitp < 0)] = np.nan
        sitp[np.where(sitp > 12)] = np.nan
    
    print 'Completed: Read PIOMAS SIT!'
    return sitp
    
#sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
#sitp = readPIOMAS(directorydatap,0.15)

lons2,lats2 = np.meshgrid(lons,lats)

sitall = np.append(sith,sitf,axis=1)

###########################################################################
###########################################################################
###########################################################################
### Calculate correlations
timesat = np.where((years >= 1979) & (years <= 2015))[0]

sitallsat = np.nanmean(sitall[:,timesat,:,:,:],axis=0)

sitsat_w = np.nanmean(sitallsat[:,0:3,:,:],axis=1)
sitsat_sp = np.nanmean(sitallsat[:,3:6,:,:],axis=1)
sitsat_su = np.nanmean(sitallsat[:,6:9,:,:],axis=1)
sitsat_f = np.nanmean(sitallsat[:,9:12,:,:],axis=1)

sitp_w = np.nanmean(sitp[:,0:3,:,:],axis=1)
sitp_sp = np.nanmean(sitp[:,3:6,:,:],axis=1)
sitp_su = np.nanmean(sitp[:,6:9,:,:],axis=1)
sitp_f = np.nanmean(sitp[:,9:12,:,:],axis=1)

def corr(sitp,sitl):
    varx = sitp
    vary = sitl
    corr = np.empty((sitp.shape[1],sitp.shape[2]))
    for i in xrange(sitp.shape[1]):
        for j in xrange(sitp.shape[2]):
            corr[i,j] = sts.stats.pearsonr(varx[:,i,j],vary[:,i,j])[0]
        
#    corr[np.where(corr == 1.)] = np.nan
    
    print 'Completed: Correlated PIOMAS and LENS SIT data!'
    return corr

### Calculate correlations
corr_w = corr(sitp_w,sitsat_w)
corr_sp = corr(sitp_sp,sitsat_sp)
corr_su = corr(sitp_su,sitsat_su)
corr_f = corr(sitp_f,sitsat_f)

### Plot figure
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

cs = m.contourf(lons2,lats2,corr_w,
                values,latlon=True)
cs1 = m.contour(lons2,lats2,corr_w,
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

cs = m.contourf(lons2,lats2,corr_sp,
                values,latlon=True)
cs1 = m.contour(lons2,lats2,corr_sp,
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

cs = m.contourf(lons2,lats2,corr_su,
                values,latlon=True)
cs1 = m.contour(lons2,lats2,corr_su,
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

cs = m.contourf(lons2,lats2,corr_f,
                values,latlon=True)
cs1 = m.contour(lons2,lats2,corr_f,
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

plt.savefig(directoryfigure + 'LENSPIOMAS_SIT_corr.png',dpi=300)