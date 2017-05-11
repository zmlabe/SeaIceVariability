"""
Script calculates and plots trends in PIOMAS SIT

Author : Zachary Labe
Date : 13 September 2016
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import scipy.stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import datetime
import iris as ir
import iris.quickplot as qplt
import read_SeaIceThick_PIOMAS as CT
import calc_PiomasArea as CA
from matplotlib.colors import Normalize
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate PIOMAS trends - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',r'Sep',
          r'Oct',r'Nov',r'Dec']
          
yearmin2 = 2005
yearmax2 = 2015
years2 = np.arange(yearmin2,yearmax2+1,1)

### Call functions
lats,lons,sitq1 = CT.readPiomas(directorydata,years,0.15)
#lats,lons,sitq2 = CT.readPiomas(directorydata,years2,0.15)
area = CA.readPiomasArea(directorydata)

#### Take monthly mean
def monRegress(sitq,months):
    slopesit = np.zeros((sitq.shape[1],sitq.shape[2],sitq.shape[3]))
    for mo in xrange(sitq.shape[1]):
        sit = sitq[:,mo,:,:]
        for i in xrange(0,sit.shape[1]):
            for j in xrange(0,sit.shape[2]):
                varyy = np.ravel(sit[:,i,j])
                varxx = np.arange(varyy.shape[0])
                mask = np.isfinite(varxx) & np.isfinite(varyy)
                
                varyymean = np.nanmean(varyy)
                if np.isfinite(varyymean):
                    slopesit[mo,i,j],intercept,r,p_value,std_err = sts.stats.linregress(varxx[mask],
                                                                      varyy[mask])
                else:
                    slopesit[mo,i,j] = np.nan  
                    
        print 'Completed: Month %s done!' % (months[mo])
    print 'Completed: Calculated regression!'
    
    slopesit = slopesit*10. # decadal trend
    return slopesit

slopesit1 = monRegress(sitq1,months)  
#slopesit2 = monRegress(sitq2,months)   
#
##slopesit1 = np.nanmean(slopesit1,axis=0)
#slopesit2 = np.nanmean(slopesit2,axis=0)
#
sittrend_w = np.nanmean(slopesit1[0:3],axis=0)
sittrend_sp = np.nanmean(slopesit1[3:6],axis=0)
sittrend_su = np.nanmean(slopesit1[6:9],axis=0)
sittrend_f = np.nanmean(slopesit1[9:12],axis=0)

def weightThick(var,area):
    """
    Area weights sit array 4d [year,month,lat,lon] into [year,month]
    """
    varq = var[:,:]
    mask = np.isfinite(varq) & np.isfinite(area)
    varmask = varq[mask]
    areamask = area[mask]
    sityr = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr
    
avetrend_w = weightThick(sittrend_w,area)
avetrend_sp = weightThick(sittrend_sp,area)
avetrend_su = weightThick(sittrend_su,area)
avetrend_f = weightThick(sittrend_f,area)

#sittrend_w = np.nanmean(slopesit2[0:3],axis=0)
#sittrend_sp = np.nanmean(slopesit2[3:6],axis=0)
#sittrend_su = np.nanmean(slopesit2[6:9],axis=0)
#sittrend_f = np.nanmean(slopesit2[9:12],axis=0)

#class MidpointNormalize(Normalize):
#    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#        self.midpoint = midpoint
#        Normalize.__init__(self, vmin, vmax, clip)
#
#    def __call__(self, value, clip=None):
#        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#        return np.ma.masked_array(np.interp(value, x, y)) 
#norm = MidpointNormalize(midpoint=0)

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

#### Define figure
#fig = plt.figure()
#for i in xrange(slopesit.shape[0]):    
#    ax = plt.subplot(3,4,i+1)
#    
#    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
#    m.drawmapboundary(fill_color='white')
#    m.drawcoastlines(color='darkgrey',linewidth=0.1)
#    parallels = np.arange(50,90,10)
#    meridians = np.arange(-180,180,30)
##    m.drawparallels(parallels,labels=[False,False,False,False],
##                    linewidth=0.0,color='k',fontsize=5)
##    m.drawmeridians(meridians,labels=[False,False,False,False],
##                    linewidth=0.0,color='k',fontsize=5)
#    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#    
#    ### Adjust maximum limits
#    values = np.arange(-0.75,0.6,0.05)  
#    
#    ### Plot filled contours    
#    cs = m.contourf(lons[:,:],lats[:,:],slopesit[i,:,:],
#                    values,latlon=True,norm=norm)
#    cs1 = m.contour(lons[:,:],lats[:,:],slopesit[i,:,:],
#                    values,latlon=True,linestyles='-',
#                    colors='k',linewidths=0.05)
#                      
#    ### Set colormap     
##    cmap = plt.cm.get_cmap('brewer_RdBu_11')      
#    cmap = plt.cm.get_cmap('bwr_r')                      
#    cs.set_cmap(cmap)
#    
#    ax.text(0.89,0.95,r'\textbf{%s}' % (months[i]),size='5',
#                horizontalalignment='center',backgroundcolor='w',
#                verticalalignment='center',bbox=dict(facecolor='w',
#                edgecolor='k',alpha=0.9),transform=ax.transAxes)
#                          
#cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='both',extendfrac=0.07,drawedges=True)
#
#cbar.set_label(r'm decade$^{-1}$')
#cbar.set_ticks(np.arange(-0.75,0.6,0.25))
#cbar.set_ticklabels(map(str,np.arange(-0.75,0.6,0.25))) 
#
#fig.suptitle(r'\textbf{PIOMAS, 1979-2015 SIT Trends}',fontsize=14)
#fig.subplots_adjust(bottom=0.15)
#fig.subplots_adjust(wspace=-0.5)
##plt.tight_layout()
#
#### Save figure
#plt.savefig(directoryfigure +'largegrid.png',dpi=500)

#fig = plt.figure()   
#ax = plt.subplot(1,2,1)
#
#m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='k',linewidth=0.3)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.35,color='k',fontsize=5)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.35,color='k',fontsize=5)
#m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#
#### Adjust maximum limits
#values = np.arange(-1,1.1,0.1) 
#
#### Plot filled contours    
#cs = m.contourf(lons[:,:],lats[:,:],slopesit1[:,:],
#                values,latlon=True,extend='both')
#cs1 = m.contour(lons[:,:],lats[:,:],slopesit1[:,:],
#                values,latlon=True,linestyles='-',
#                colors='k',linewidths=0.05)
#                  
#### Set colormap         
#cmap = plt.cm.get_cmap('bwr_r')                      
#cs.set_cmap(cmap)
#
#ax = plt.subplot(1,2,2)
#
#m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='k',linewidth=0.3)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.35,color='k',fontsize=5)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.35,color='k',fontsize=5)
#m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#
#### Adjust maximum limits
#values = np.arange(-1,1.1,.1)
#
#### Plot filled contours    
#cs = m.contourf(lons[:,:],lats[:,:],slopesit2[:,:],
#                values,latlon=True,extend='both')
#cs1 = m.contour(lons[:,:],lats[:,:],slopesit2[:,:],
#                values,latlon=True,linestyles='-',
#                colors='k',linewidths=0.05)
#                  
#### Set colormap     
##    cmap = plt.cm.get_cmap('brewer_RdBu_11')      
#cmap = plt.cm.get_cmap('seismic_r')                      
#cs.set_cmap(cmap)
#                          
#cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='both',extendfrac=0.07,drawedges=True)
#
#cbar.set_label(r'SIT( m decade$^{-1}$ )')
#cbar.set_ticks(np.arange(-1,1.1,0.5))
#cbar.set_ticklabels(map(str,np.arange(-1,1.1,0.5))) 
#
#plt.text(0.82,26,'2000-2015',fontsize=23)
#plt.text(-0.56,26,'1979-1999',fontsize=23)
#
##fig.subplots_adjust(bottom=0.15)
##fig.subplots_adjust(wspace=-0.5)
#plt.tight_layout()
#
#### Save figure
#plt.savefig(directoryfigure +'largegrid_sitsubs.png',dpi=500)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
### Plot Composites
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

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
barlim = np.arange(-1,2,1)
values = np.arange(-1,1.1,0.1)

cs = m.contourf(lons,lats,var,
                values,latlon=True,extend='both')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BlueDarkOrange18')         
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
                values,latlon=True,extend='both')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BlueDarkOrange18')         
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
                values,latlon=True,extend='both')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BlueDarkOrange18')         
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
                values,latlon=True,extend='both')
cs1 = m.contour(lons,lats,var,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-',latlon=True)
        
cmap = ncm.cmap('BlueDarkOrange18')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='Both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{SIT( m decade$^{-1}$ )}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'sit_trends_7915.png',dpi=300)

print 'Completed: Script done!'