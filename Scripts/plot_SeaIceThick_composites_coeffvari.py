"""
Scripts plots sit from LENS future
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 16 November 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import datetime
import read_SeaIceThick_LENS as lens
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
from netCDF4 import Dataset

### Define directories
directorydatal = '/home/zlabe/Surtsey3/'
directorydatap = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----LENS Historical Mean Sea Ice Thickness - %s----' % titletime 

### Alott time series
yearmin = 1920
yearmax = 2080
years = np.arange(yearmin,yearmax+1,1)
years2 = np.arange(2006,2080+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,39,1)) + map(str,np.arange(101,106,1))

def weightThick(var,lats,types):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    """
    
    if types == 'lens':
        sityr = np.empty((var.shape[0],var.shape[1],var.shape[2]))
        for ens in xrange(var.shape[0]):
            for i in xrange(var.shape[1]):
                for j in xrange(var.shape[2]):
                    varq = var[ens,i,j,:,:]
                    mask = np.isfinite(varq) & np.isfinite(lats)
                    varmask = varq[mask]
                    areamask = np.cos(np.deg2rad(lats[mask]))
                    sityr[ens,i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
            
            print 'Completed: Weighting per ensemble #%s!' % ensemble[ens]
    
    elif types == 'piomas':
        sityr = np.empty((var.shape[0],var.shape[1]))
        for i in xrange(var.shape[0]):
            for j in xrange(var.shape[1]):
                varq = var[i,j,:,:]
                mask = np.isfinite(varq) & np.isfinite(lats)
                varmask = varq[mask]
                areamask = np.cos(np.deg2rad(lats[mask]))
                sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr

### Call functions   
#sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
#lons2,lats2 = np.meshgrid(lons,lats)
#  
#sitaveh = weightThick(sith,lats2,'lens')
#sitavef = weightThick(sitf,lats2,'lens')

yearp1 = np.where((years >= 1980) & (years <= 1997))[0]
yearp2 = np.where((years >= 1998) & (years <= 2015))[0]
yearqh1 = np.where((years >= 1920) & (years <= 1962))[0]
yearqh2 = np.where((years >= 1963) & (years <= 2005))[0]
yearqf1 = np.where((years2 >= 2006) & (years2 <= 2042))[0]
yearqf2 = np.where((years2 >= 2043) & (years2 <= 2080))[0]

#### September 
sith_mo2 = sith[:,:,8,:,:]
sitf_mo2 = sitf[:,:,8,:,:]
sitall_mo2 = np.append(sith[:,:,8,:,:],sitf[:,:,8,:,:],axis=1)

sith1 = np.nanmean(sith_mo2[:,yearqh1],axis=1)
sith2 = np.nanmean(sith_mo2[:,yearqh2],axis=1)
sith3 = np.nanmean(sitall_mo2[:,yearp1],axis=1)

sitf1 = np.nanmean(sitf_mo2[:,yearqf1,:,:],axis=1)
sitf2 = np.nanmean(sitf_mo2[:,yearqf2,:,:],axis=1)
sitf3 = np.nanmean(sitall_mo2[:,yearp2,:,:],axis=1)

### Max/min ensembles
sith1diff = (np.nanstd(sith1,axis=0)/np.nanmean(sith1,axis=0)) * 100.
sith2diff = (np.nanstd(sith2,axis=0)/np.nanmean(sith2,axis=0)) * 100.
sith3diff = (np.nanstd(sith3,axis=0)/np.nanmean(sith3,axis=0)) * 100.

sitf1diff = (np.nanstd(sitf1,axis=0)/np.nanmean(sitf1,axis=0)) * 100.
sitf2diff = (np.nanstd(sitf2,axis=0)/np.nanmean(sitf2,axis=0)) * 100.
sitf3diff = (np.nanstd(sitf3,axis=0)/np.nanmean(sitf3,axis=0)) * 100.

#sith1diff = np.nanstd(sith1,axis=0)
#sith2diff = np.nanstd(sith2,axis=0)
#sith3diff = np.nanstd(sith3,axis=0)
#
#sitf1diff = np.nanstd(sitf1,axis=0)
#sitf2diff = np.nanstd(sitf2,axis=0)
#sitf3diff = np.nanstd(sitf3,axis=0)

#sith1diff = np.nanmax(sith1,axis=0) - np.nanmin(sith1,axis=0)
#sith2diff = np.nanmax(sith2,axis=0) - np.nanmin(sith2,axis=0)
#sith3diff = np.nanmax(sith3,axis=0) - np.nanmin(sith3,axis=0)
#
#sitf1diff = np.nanmax(sitf1,axis=0) - np.nanmin(sitf1,axis=0)
#sitf2diff = np.nanmax(sitf2,axis=0) - np.nanmin(sitf2,axis=0)
#sitf3diff = np.nanmax(sitf3,axis=0) - np.nanmin(sitf3,axis=0)

composites = [sith1diff,sith2diff,sith3diff,sitf1diff,sitf2diff,sitf3diff]

### Create subplots
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

fig = plt.figure()
    
for i in xrange(len(composites)):
    ax = plt.subplot(2,3,i+1)
    
    ### Select variable
    var = composites[i]
    
    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
                resolution='l',round =True)
                
    var, lons_cyclic = addcyclic(var, lons)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
    x, y = m(lon2d, lat2d)      
      
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.2)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0.35,color='k',fontsize=1)
#    m.drawmeridians(meridians,labels=[False,False,False,False],
#                    linewidth=0.35,color='k',fontsize=1)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(0,101,10)  
#    values = np.arange(0,1.1,0.1) 
#    values = np.arange(-3,3.1,0.25)  
    
    ### Plot filled contours    
    cs = m.contourf(x,y,var,
                    values,extend='both')
    cs1 = m.contour(x,y,var,
                    values,linewidths=0.2,colors='darkgrey',
                    linestyles='-')
                    
    ### Set colormap  
    cmap = plt.cm.get_cmap('cubehelix_r')         
    cs.set_cmap(cmap) 
    
cbar_ax = fig.add_axes([0.313,0.13,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)  

perc = r'$\bf{\%}$'                    
cbar.set_ticks(np.arange(0,101,25))
cbar.set_ticklabels(map(str,np.arange(0,101,25)))   
#cbar.set_ticks(np.arange(0,1.1,0.5))
#cbar.set_ticklabels(map(str,np.arange(0,1.1,0.5)))  
#cbar.set_ticks(np.arange(-3,4,1))
#cbar.set_ticklabels(map(str,np.arange(-3,4,1)))  
cbar.set_label(r'\textbf{CV (%s)}' % perc)
cbar.ax.tick_params(axis='x', size=.01)

plt.subplots_adjust(wspace=-0.28)
plt.subplots_adjust(hspace=0.15)
plt.subplots_adjust(bottom=0.2)
plt.subplots_adjust(top=0.87)
            
plt.annotate(r'\textbf{LENS}', xy=(0, 0), xytext=(0.35, 0.915),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0)
plt.annotate(r'\textbf{LENS}', xy=(0, 0), xytext=(0.680, 0.915),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0)
plt.annotate(r'\textbf{PIOMAS}', xy=(0, 0), xytext=(0.773, 0.892),
            xycoords='figure fraction',fontsize=7,color='darkgrey',
            rotation=0,ha='center')            
            
plt.annotate(r'\textbf{1980-1997}', xy=(0, 0), xytext=(0.79, 0.863),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)                  
plt.annotate(r'\textbf{1998-2015}', xy=(0, 0), xytext=(0.79, 0.504),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)

plt.annotate(r'\textbf{1920-1962}', xy=(0, 0), xytext=(0.24, 0.875),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)                               
plt.annotate(r'\textbf{1963-2005}', xy=(0, 0), xytext=(0.468, 0.875),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)
plt.annotate(r'\textbf{2006-2042}', xy=(0, 0), xytext=(0.24, 0.517),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)                     
plt.annotate(r'\textbf{2043-2080}', xy=(0, 0), xytext=(0.468, 0.517),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0) 
    
### Save figure
plt.savefig(directoryfigure +'sit_rcp_composites_coeffvari.png',dpi=500)