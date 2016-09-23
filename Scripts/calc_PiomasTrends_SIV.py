"""
Script calculates and plots trends in PIOMAS SIV

Author : Zachary Labe
Date : 13 September 2016
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import scipy.stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
import iris as ir
import iris.quickplot as qplt
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
import calc_SIV as CV
from matplotlib.colors import Normalize

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
yearmax = 1999
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',r'Sep',
          r'Oct',r'Nov',r'Dec']
          
yearmin2 = 2000
yearmax2 = 2015
years2 = np.arange(yearmin2,yearmax2+1,1)

### Call functions
#lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
#lats,lons,sic = CC.readPiomas(directorydata,years,0.01)
#area = CA.readPiomasArea(directorydata)
#sivq = CV.sivGrid(sit,sic,area,True)

#lats,lons,sit2 = CT.readPiomas(directorydata,years2,0.15)
#lats,lons,sic2 = CC.readPiomas(directorydata,years2,0.01)
#sivq2 = CV.sivGrid(sit2,sic2,area,True)

### Take monthly mean
def monRegress(sivq,months):
    slopesiv = np.zeros((sivq.shape[1],sivq.shape[2],sivq.shape[3]))
    for mo in xrange(sivq.shape[1]):
        siv = sivq[:,mo,:,:]
        for i in xrange(0,siv.shape[1]):
            for j in xrange(0,siv.shape[2]):
                varyy = np.ravel(siv[:,i,j])
                varxx = np.arange(varyy.shape[0])
                mask = np.isfinite(varxx) & np.isfinite(varyy)
                
                varyymean = np.nanmean(varyy)
                if np.isfinite(varyymean):
                    slopesiv[mo,i,j],intercept,r,p_value,std_err = sts.stats.linregress(varxx[mask],
                                                                      varyy[mask])
                else:
                    slopesiv[mo,i,j] = np.nan  
        print 'Completed: Month %s done!' % (months[mo])
    print 'Completed: Calculated regression!'
    
    slopesiv = slopesiv*10. # decadal trend
    return slopesiv

#slopesiv1 = monRegress(sivq,months)  
#slopesiv2 = monRegress(sivq2,months)  

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

#### Define figure
#fig = plt.figure()
#for i in xrange(slopesiv.shape[0]):    
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
#    values = np.arange(-1000,1001,50)  
#    
#    ### Plot filled contours    
#    cs = m.contourf(lons[:,:],lats[:,:],slopesiv[i,:,:],
#                    values,latlon=True,extend='both')
#    cs1 = m.contour(lons[:,:],lats[:,:],slopesiv[i,:,:],
#                    values,latlon=True,linestyles='-',
#                    colors='k',linewidths=0.05)
#                      
#    ### Set colormap     
##    cmap = plt.cm.get_cmap('brewer_RdBu_11')      
#    cmap = plt.cm.get_cmap('spectral')                      
#    cs.set_cmap(cmap)
#    
#    ax.text(0.91,0.97,r'\textbf{%s}' % (months[i]),size='8',
#                horizontalalignment='center',
#                verticalalignment='center',transform=ax.transAxes)
#                          
#cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='both',extendfrac=0.07,drawedges=True)
#
#cbar.set_label(r'km$^{3}$ decade$^{-1}$')
#cbar.set_ticks(np.arange(-1000,1001,500))
#cbar.set_ticklabels(map(str,np.arange(-1000,1001,500))) 
#
#fig.suptitle(r'\textbf{PIOMAS, 1979-2015 SIV Trends}',fontsize=14)
#fig.subplots_adjust(bottom=0.15)
#fig.subplots_adjust(wspace=-0.5)
##plt.tight_layout()
#
#### Save figure
#plt.savefig(directoryfigure +'largegrid_vol.png',dpi=500)

slopesiv1 = np.nanmean(slopesiv1,axis=0)
slopesiv2 = np.nanmean(slopesiv2,axis=0)

fig = plt.figure()   
ax = plt.subplot(1,2,1)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.35,color='k',fontsize=5)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.35,color='k',fontsize=5)
m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')

### Adjust maximum limits
values = np.arange(-1000,1001,50)  

### Plot filled contours    
cs = m.contourf(lons[:,:],lats[:,:],slopesiv1[:,:],
                values,latlon=True,extend='both')
cs1 = m.contour(lons[:,:],lats[:,:],slopesiv1[:,:],
                values,latlon=True,linestyles='-',
                colors='k',linewidths=0.05)
                  
### Set colormap         
cmap = plt.cm.get_cmap('seismic_r')                      
cs.set_cmap(cmap)

ax = plt.subplot(1,2,2)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],
                linewidth=0.35,color='k',fontsize=5)
m.drawmeridians(meridians,labels=[False,False,False,False],
                linewidth=0.35,color='k',fontsize=5)
m.drawlsmask(land_color='dimgrey',ocean_color='mintcream')

### Adjust maximum limits
values = np.arange(-1000,1001,50)  

### Plot filled contours    
cs = m.contourf(lons[:,:],lats[:,:],slopesiv2[:,:],
                values,latlon=True,extend='both')
cs1 = m.contour(lons[:,:],lats[:,:],slopesiv2[:,:],
                values,latlon=True,linestyles='-',
                colors='k',linewidths=0.05)
                  
### Set colormap     
#    cmap = plt.cm.get_cmap('brewer_RdBu_11')      
cmap = plt.cm.get_cmap('seismic_r')                      
cs.set_cmap(cmap)
                          
cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'SIV( km$^{3}$ decade$^{-1}$ )')
cbar.set_ticks(np.arange(-1000,1001,500))
cbar.set_ticklabels(map(str,np.arange(-1000,1001,500))) 

plt.text(0.82,26,'2000-2015',fontsize=23)
plt.text(-0.56,26,'1979-1999',fontsize=23)

#fig.subplots_adjust(bottom=0.15)
#fig.subplots_adjust(wspace=-0.5)
plt.tight_layout()

### Save figure
plt.savefig(directoryfigure +'largegrid_volsubs.png',dpi=500)

print 'Completed: Script done!'