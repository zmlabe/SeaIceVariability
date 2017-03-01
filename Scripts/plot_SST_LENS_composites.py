"""
Scripts plots SST for CESM Large ensemble over historical and future
periods

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 8 December 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_var_LENS as LV
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey3/CESM_large_ensemble/' 
directoryfigure = '/home/zlabe/Desktop/LENS/SST/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot LENS AD - %s----' % titletime 

ensembles = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 1980
year2 = 2015
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

yearslens = np.arange(1920,2080+1,1)          
yearsclimo = np.arange(1981,2010+1,1)
          
### Read in functions
sstk,lats,lons = LV.readLENSEnsemble(directorydata,'SST') 

### Pick years
yearq = np.where((yearslens >= 2006) & (yearslens <= 2080))[0]
sstk = sstk[:,yearq,:,:,:]

### Mesh lat/lon
#lons,lats = np.meshgrid(lons,lats)

### Remove missing data (<0 Kelvin!!)
sstk[np.where(sstk < 0)] = np.nan

## Convert to celsius 
sst = sstk - 273.15 

### Mask freezing
sst[np.where(sst <= -1.8)] = -1.8 # freezing point SST-ice

### Calculate seasons
sst_w = np.nanmean(sst[:,:,0:3,:,:],axis=2)
sst_sp = np.nanmean(sst[:,:,3:6,:,:],axis=2)
sst_su = np.nanmean(sst[:,:,6:9,:,:],axis =2)
sst_f = np.nanmean(sst[:,:,9:12,:,:],axis=2)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Calculate decadal trend

### Calculate seasonal trend
def Regress(sstq):
    slopesst = np.zeros((sstq.shape[1],sstq.shape[2]))
    for i in xrange(0,sstq.shape[1]):
        for j in xrange(0,sstq.shape[2]):
            varyy = np.ravel(sstq[:,i,j])
            varxx = np.arange(varyy.shape[0])
            mask = np.isfinite(varxx) & np.isfinite(varyy)
            
            varyymean = np.nanmean(varyy)
            if np.isfinite(varyymean):
                slopesst[i,j],intercept,r,p_value,std_err = sts.stats.linregress(varxx[mask],
                                                                  varyy[mask])
            else:
                slopesst[i,j] = np.nan  
    print 'Completed: Calculated regression!'
    
    slopesst = slopesst*10. # decadal trend
    return slopesst

slopesst_wn = []
slopesst_spn = []
slopesst_sun = []
slopesst_fn = []
for i in xrange(sst.shape[0]):
    slopesst_wq = Regress(sst_w[i]) 
    slopesst_spq = Regress(sst_sp[i]) 
    slopesst_suq = Regress(sst_su[i]) 
    slopesst_fq = Regress(sst_f[i]) 
    
    slopesst_wn.append(slopesst_wq)
    slopesst_spn.append(slopesst_spq)
    slopesst_sun.append(slopesst_suq)
    slopesst_fn.append(slopesst_fq)
    
    print 'Completed: Regressing Ensemble #%s!' % ensembles[i]

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Calculate ensemble means
    
slopesst_w = np.nanmean(np.asarray(slopesst_wn),axis=0)
slopesst_sp = np.nanmean(np.asarray(slopesst_spn),axis=0)
slopesst_su = np.nanmean(np.asarray(slopesst_sun),axis=0)
slopesst_f = np.nanmean(np.asarray(slopesst_fn),axis=0)

#### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

### Append seasons into list
sstrend_seasons = [slopesst_w,slopesst_sp,slopesst_su,slopesst_f]

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Plot Composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=50,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(slopesst_w, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

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

cs = m.contourf(x,y,var,
                values)
cs1 = m.contour(x,y,var,
                values,linewidths=0.2,colors='k',
                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JFM}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=50,lon_0=270,
            resolution='l',round =True)

var, lons_cyclic = addcyclic(slopesst_sp, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,
                values)
cs1 = m.contour(x,y,var,
                values,linewidths=0.2,colors='k',
                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(223)

m = Basemap(projection='npstere',boundinglat=50,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(slopesst_su, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
                    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,
                values)
cs1 = m.contour(x,y,var,
                values,linewidths=0.2,colors='k',
                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(224)

m = Basemap(projection='npstere',boundinglat=50,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(slopesst_f, lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,
                values,extend='both')
cs1 = m.contour(x,y,var,
                values,linewidths=0.2,colors='k',
                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='Both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{Decadal Trend SST($^{\circ}$C)}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'sst_seasonal_trends_rcp85.png',dpi=300)