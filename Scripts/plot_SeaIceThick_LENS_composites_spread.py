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
yearmin = 2006
yearmax = 2080
years = np.arange(yearmin,yearmax+1,1)
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

### Calculate sit_max and sit_min
sitaveh1 = np.nanmean(sitaveh[:,11:36,8],axis=1)
sitaveh2 = np.nanmean(sitaveh[:,36:61,8],axis=1)
sitaveh3 = np.nanmean(sitaveh[:,61:86,8],axis=1)

sitavef1 = np.nanmean(sitavef[:,0:25,8],axis=1)
sitavef2 = np.nanmean(sitavef[:,25:50,8],axis=1)
sitavef3 = np.nanmean(sitavef[:,50:75,8],axis=1)

maxh1 = np.where(sitaveh1 == np.nanmax(sitaveh1))[0]
maxh2 = np.where(sitaveh2 == np.nanmax(sitaveh2))[0]
maxh3 = np.where(sitaveh3 == np.nanmax(sitaveh3))[0]

minh1 = np.where(sitaveh1 == np.nanmin(sitaveh1))[0]
minh2 = np.where(sitaveh2 == np.nanmin(sitaveh2))[0]
minh3 = np.where(sitaveh3 == np.nanmin(sitaveh3))[0]

maxf1 = np.where(sitavef1 == np.nanmax(sitavef1))[0]
maxf2 = np.where(sitavef2 == np.nanmax(sitavef2))[0]
maxf3 = np.where(sitavef3 == np.nanmax(sitavef3))[0]

minf1 = np.where(sitavef1 == np.nanmin(sitavef1))[0]
minf2 = np.where(sitavef2 == np.nanmin(sitavef2))[0]
minf3 = np.where(sitavef3 == np.nanmin(sitavef3))[0]

#### September 
sith_mo = sith[:,:,8,:,:]

sith1 = np.nanmean(sith_mo[:,11:36],axis=1)
sith2 = np.nanmean(sith_mo[:,36:61],axis=1)
sith3 = np.nanmean(sith_mo[:,61:86],axis=1)

sitf_mo = sitf[:,:,8,:,:]

sitf1 = np.nanmean(sitf_mo[:,0:25],axis=1)
sitf2 = np.nanmean(sitf_mo[:,25:50],axis=1)
sitf3 = np.nanmean(sitf_mo[:,50:75],axis=1)

### Max/min ensembles
#sith1diff = np.squeeze(sith1[maxh1] - sith1[minh1])
#sith2diff = np.squeeze(sith2[maxh2] - sith2[minh2])
#sith3diff = np.squeeze(sith3[maxh3] - sith3[minh3])
#
#sitf1diff = np.squeeze(sitf1[maxf1] - sitf1[minf1])
#sitf2diff = np.squeeze(sitf2[maxf2] - sitf2[minf2])
#sitf3diff = np.squeeze(sitf3[maxf3] - sitf3[minf3])

### Percentiles 90/10
#sortedh1 = np.sort(sitaveh1)
#sortedh2 = np.sort(sitaveh2)
#sortedh3 = np.sort(sitaveh3)
#
#sortedf1 = np.sort(sitavef1)
#sortedf2 = np.sort(sitavef2)
#sortedf3 = np.sort(sitavef3)
#
#mxh1 = np.where(sitaveh1 == sortedh1[-4])
#mxh2 = np.where(sitaveh2 == sortedh2[-4])
#mxh3 = np.where(sitaveh3 == sortedh3[-4])
#
#mnh1 = np.where(sitaveh1 == sortedh1[3])
#mnh2 = np.where(sitaveh2 == sortedh2[3])
#mnh3 = np.where(sitaveh3 == sortedh3[3])
#
#mxf1 = np.where(sitavef1 == sortedf1[-4])
#mxf2 = np.where(sitavef2 == sortedf2[-4])
#mxf3 = np.where(sitavef3 == sortedf3[-4])
#
#mnf1 = np.where(sitavef1 == sortedf1[3])
#mnf2 = np.where(sitavef2 == sortedf2[3])
#mnf3 = np.where(sitavef3 == sortedf3[3])

### Percentiles 75/25
#sortedh1 = np.sort(sitaveh1)
#sortedh2 = np.sort(sitaveh2)
#sortedh3 = np.sort(sitaveh3)
#
#sortedf1 = np.sort(sitavef1)
#sortedf2 = np.sort(sitavef2)
#sortedf3 = np.sort(sitavef3)
#
#mxh1 = np.where(sitaveh1 == sortedh1[-10])
#mxh2 = np.where(sitaveh2 == sortedh2[-10])
#mxh3 = np.where(sitaveh3 == sortedh3[-10])
#
#mnh1 = np.where(sitaveh1 == sortedh1[9])
#mnh2 = np.where(sitaveh2 == sortedh2[9])
#mnh3 = np.where(sitaveh3 == sortedh3[9])
#
#mxf1 = np.where(sitavef1 == sortedf1[-10])
#mxf2 = np.where(sitavef2 == sortedf2[-10])
#mxf3 = np.where(sitavef3 == sortedf3[-10])
#
#mnf1 = np.where(sitavef1 == sortedf1[9])
#mnf2 = np.where(sitavef2 == sortedf2[9])
#mnf3 = np.where(sitavef3 == sortedf3[9])
#
#sith1diff = np.squeeze(sith1[mxh1] - sith1[mnh1])
#sith2diff = np.squeeze(sith2[mxh2] - sith2[mnh2])
#sith3diff = np.squeeze(sith3[mxh3] - sith3[mnh3])
#
#sitf1diff = np.squeeze(sitf1[mxf1] - sitf1[mnf1])
#sitf2diff = np.squeeze(sitf2[mxf2] - sitf2[mnf2])
#sitf3diff = np.squeeze(sitf3[mxf3] - sitf3[mnf3])

#### Percentiles 90/10 grid cell method
#p90h_1 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10h_1 = np.empty((sith1.shape[1],sith1.shape[2]))
#p90h_2 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10h_2 = np.empty((sith1.shape[1],sith1.shape[2]))
#p90h_3 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10h_3 = np.empty((sith1.shape[1],sith1.shape[2]))
#
#p90f_1 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10f_1 = np.empty((sith1.shape[1],sith1.shape[2]))
#p90f_2 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10f_2 = np.empty((sith1.shape[1],sith1.shape[2]))
#p90f_3 = np.empty((sith1.shape[1],sith1.shape[2]))
#p10f_3 = np.empty((sith1.shape[1],sith1.shape[2]))
#
#for i in xrange(sith1.shape[1]):
#    for j in xrange(sith1.shape[2]):
#        p90h_1[i,j] = np.nanpercentile(sith1[:,i,j],90)
#        p90h_2[i,j] = np.nanpercentile(sith2[:,i,j],90)
#        p90h_3[i,j] = np.nanpercentile(sith3[:,i,j],90)
#        
#        p10h_1[i,j] = np.nanpercentile(sith1[:,i,j],10)
#        p10h_2[i,j] = np.nanpercentile(sith2[:,i,j],10)
#        p10h_3[i,j] = np.nanpercentile(sith3[:,i,j],10)
#        
#        p90f_1[i,j] = np.nanpercentile(sitf1[:,i,j],90)
#        p90f_2[i,j] = np.nanpercentile(sitf2[:,i,j],90)
#        p90f_3[i,j] = np.nanpercentile(sitf3[:,i,j],90)
#        
#        p10f_1[i,j] = np.nanpercentile(sitf1[:,i,j],10)
#        p10f_2[i,j] = np.nanpercentile(sitf2[:,i,j],10)
#        p10f_3[i,j] = np.nanpercentile(sitf3[:,i,j],10)
#        
#sith1diff = p90h_1 - p10h_1
#sith2diff = p90h_2 - p10h_2
#sith3diff = p90h_3 - p10h_3
#
#sitf1diff = p90f_1 - p10f_1
#sitf2diff = p90f_2 - p10f_2
#sitf3diff = p90f_3 - p10f_3

### Max/min grid cell method
p90h_1 = np.empty((sith1.shape[1],sith1.shape[2]))
p10h_1 = np.empty((sith1.shape[1],sith1.shape[2]))
p90h_2 = np.empty((sith1.shape[1],sith1.shape[2]))
p10h_2 = np.empty((sith1.shape[1],sith1.shape[2]))
p90h_3 = np.empty((sith1.shape[1],sith1.shape[2]))
p10h_3 = np.empty((sith1.shape[1],sith1.shape[2]))

p90f_1 = np.empty((sith1.shape[1],sith1.shape[2]))
p10f_1 = np.empty((sith1.shape[1],sith1.shape[2]))
p90f_2 = np.empty((sith1.shape[1],sith1.shape[2]))
p10f_2 = np.empty((sith1.shape[1],sith1.shape[2]))
p90f_3 = np.empty((sith1.shape[1],sith1.shape[2]))
p10f_3 = np.empty((sith1.shape[1],sith1.shape[2]))

for i in xrange(sith1.shape[1]):
    for j in xrange(sith1.shape[2]):
        p90h_1[i,j] = np.nanmax(sith1[:,i,j])
        p90h_2[i,j] = np.nanmax(sith2[:,i,j])
        p90h_3[i,j] = np.nanmax(sith3[:,i,j])
        
        p10h_1[i,j] = np.nanmin(sith1[:,i,j])
        p10h_2[i,j] = np.nanmin(sith2[:,i,j])
        p10h_3[i,j] = np.nanmin(sith3[:,i,j])
        
        p90f_1[i,j] = np.nanmax(sitf1[:,i,j])
        p90f_2[i,j] = np.nanmax(sitf2[:,i,j])
        p90f_3[i,j] = np.nanmax(sitf3[:,i,j])
        
        p10f_1[i,j] = np.nanmin(sitf1[:,i,j])
        p10f_2[i,j] = np.nanmin(sitf2[:,i,j])
        p10f_3[i,j] = np.nanmin(sitf3[:,i,j])
        
sith1diff = p90h_1 - p10h_1
sith2diff = p90h_2 - p10h_2
sith3diff = p90h_3 - p10h_3

sitf1diff = p90f_1 - p10f_1
sitf2diff = p90f_2 - p10f_2
sitf3diff = p90f_3 - p10f_3

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
    values = np.arange(-3,3.1,0.25)  
    
    ### Plot filled contours    
    cs = m.contourf(x,y,var,
                    values,extend='both')
    cs1 = m.contour(x,y,var,
                    values,linewidths=0.2,colors='darkgrey',
                    linestyles='-')
                    
    ### Set colormap  
    cmap = ncm.cmap('temp_19lev')        
    cs.set_cmap(cmap) 
    
cbar_ax = fig.add_axes([0.313,0.13,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)  
                    
cbar.set_ticks(np.arange(-3,4,1))
cbar.set_ticklabels(map(str,np.arange(-3,4,1)))    
cbar.set_label(r'\textbf{Difference (meters)}')

plt.subplots_adjust(hspace=0.05)
plt.subplots_adjust(wspace=-0)
plt.subplots_adjust(bottom=0.2)

plt.annotate(r'\textbf{Historical}', xy=(0, 0), xytext=(0.065, 0.845),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=90)
plt.annotate(r'\textbf{RCP8.5}', xy=(0, 0), xytext=(0.065, 0.442),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=90)
            
plt.annotate(r'\textbf{1}', xy=(0, 0), xytext=(0.244, 0.915),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0)
plt.annotate(r'\textbf{2}', xy=(0, 0), xytext=(0.5, 0.915),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0)
plt.annotate(r'\textbf{3}', xy=(0, 0), xytext=(0.756, 0.915),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0)
    
### Save figure
plt.savefig(directoryfigure +'sit_rcp_composites_spread_maxmingridcell.png',dpi=500)