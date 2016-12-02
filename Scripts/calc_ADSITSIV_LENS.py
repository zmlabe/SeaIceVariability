"""
Scripts calculates correlations between SIT and AD for LENS
 
Notes
-----
    Source : LENS
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 29 November 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from scipy import signal
import read_SeaIceThick_LENS as lens
import nclcmaps as ncm

from mpl_toolkits.basemap import Basemap

### Define directories
directorydata = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/' 
directorydata3 = '/home/zlabe/Surtsey3/' 
directoryfigure = '/home/zlabe/Desktop/LENS/AD/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot NCEP - %s----' % titletime 

### Alott time series
year1 = 1920
year2 = 2080
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in data
seasons = ['JFM','AMJ','JAS','OND'] 
ens = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1)) 

ad = np.empty((4,len(ens),len(years)))
for i in xrange(len(seasons)):
    filename1 = 'AD_%s_%s%s_LENS.txt' % (seasons[i],year1,year2)
    ad[i,:,:] = np.genfromtxt(directorydata + filename1)
    
    print 'Reading data file %s!' % seasons[i]
    
sith,lats,lons = lens.readLENSEnsemble(directorydata3,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydata3,0.15,'rcp85')

##########################################################################
##########################################################################
##########################################################################
# Calculate seasons 1920-2080
sitnew = np.append(sith,sitf,axis=1)
    
sit_w = np.nanmean(sitnew[:,:,0:3,:,:],axis=2)
ad_w = ad[0]
sit_sp = np.nanmean(sitnew[:,:,3:6,:,:],axis=2)
ad_sp = ad[1]
sit_su = np.nanmean(sitnew[:,:,6:9,:,:],axis=2)
ad_su = ad[2]
sit_f = np.nanmean(sitnew[:,:,9:12,:,:],axis=2)
ad_f = ad[3]

print 'Completed: Added SIT and calculated seasons!'

############################################################################
############################################################################
############################################################################

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
        y_detrend[i,:,:] = y[i,:,:] - (slopes*x[i] + intercepts)
     
    print 'Completed: Detrended SIT data!' 
    return y_detrend

sit_w_dtn = [] 
sit_sp_dtn = [] 
sit_su_dtn = [] 
sit_f_dtn = [] 
for i in xrange(len(ens)):
    sit_w_dtq = deTrend(sit_w[i])
    sit_sp_dtq = deTrend(sit_sp[i])
    sit_su_dtq = deTrend(sit_su[i])
    sit_f_dtq = deTrend(sit_f[i])
    
    sit_w_dtn.append(sit_w_dtq)
    sit_sp_dtn.append(sit_sp_dtq)
    sit_su_dtn.append(sit_su_dtq)
    sit_f_dtn.append(sit_f_dtq)
    
    print 'Completed: Detrend ensemble #%s!' % ens[i]
    
sit_w_dt = np.asarray(sit_w_dtn)
sit_sp_dt = np.asarray(sit_sp_dtn)
sit_su_dt = np.asarray(sit_su_dtn)
sit_f_dt = np.asarray(sit_f_dtn)

#############################################################################
#############################################################################
#############################################################################

def corr(sivgrq,adq):
    varx = sivgrq
    vary = adq
    corr = np.empty((sivgrq.shape[1],sivgrq.shape[2]))
    for i in xrange(sivgrq.shape[1]):
        for j in xrange(sivgrq.shape[2]):
            corr[i,j] = sts.stats.pearsonr(varx[:,i,j],vary)[0]
        
    corr[np.where(corr == 1.)] = np.nan
    
    print 'Completed: Correlated SIV and AD data!'
    return corr

### detrend 
corr_wn = []
corr_spn = []
corr_sun = []
corr_fn = []

for i in xrange(len(ens)):   
    corr_wq = corr(sit_w_dt[i],ad_w[i])
    corr_spq = corr(sit_sp_dt[i],ad_sp[i])
    corr_suq = corr(sit_su_dt[i],ad_su[i])
    corr_fq = corr(sit_f_dt[i],ad_f[i])
    
    corr_wn.append(corr_wq)
    corr_spn.append(corr_spq)
    corr_sun.append(corr_suq)
    corr_fn.append(corr_fq)
    
corr_w = np.asarray(corr_wn)
corr_sp = np.asarray(corr_spn)
corr_su = np.asarray(corr_sun)
corr_f = np.asarray(corr_fn)
    
corr_w = np.nanmean(np.asarray(corr_wn),axis=0)
corr_sp = np.nanmean(np.asarray(corr_spn),axis=0)
corr_su = np.nanmean(np.asarray(corr_sun),axis=0)
corr_f = np.nanmean(np.asarray(corr_fn),axis=0)    

###########################################################################
###########################################################################
###########################################################################
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

lons,lats = np.meshgrid(lons,lats)

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-1,1.1,.5)
values = np.arange(-1,1.1,0.1)

cs = m.contourf(lons,lats,corr_w,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_w,
                values,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BrownBlue12')        
cs.set_cmap(cmap) 
ax.annotate(r'\textbf{JFM}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22)

###########################################################################
###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_sp,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_sp,
                values,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BrownBlue12')        
cs.set_cmap(cmap) 

ax.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22)

###########################################################################
###########################################################################

ax = plt.subplot(223)

m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_su,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_su,
                values,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BrownBlue12')        
cs.set_cmap(cmap) 
ax.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22)

###########################################################################
###########################################################################

ax = plt.subplot(224)

m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(lons,lats,corr_f,
                values,latlon=True)
cs1 = m.contour(lons,lats,corr_f,
                values,linewidths=0.2,colors='k',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BrownBlue12')        
cs.set_cmap(cmap) 

ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22)

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

plt.savefig(directoryfigure + 'lens_adsit_dt',dpi=300)
    
    