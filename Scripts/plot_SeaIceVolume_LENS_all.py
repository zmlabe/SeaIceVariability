"""
Scripts calculates SIT from LENS
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 20 October 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import datetime
import read_SeaIceThick_LENS as lens
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
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
print '\n' '----LENS Mean Sea Ice Volume - %s----' % titletime 

### Alott time series
yearmin = 1920
yearmax = 2100
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ens = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,29,1)) 

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
    
def areaGrid(lat1,lon1):
    import math as m
    
    lon2,lat2 = np.meshgrid(lon1,lat1)
    
    lat2[:] = 110.54 #km
    latsize = lat2
    
    lonsize = np.empty((lon2.shape))
    for i in xrange(len(lats)):
        lonsize[i,:] = abs(111.320*np.cos(m.radians(lat1[i])))
        
    area = latsize*lonsize
    
    print 'Completed: Calculated area of grid!'
    return area,lonsize,latsize
    
### Calculate volume
def volumeYear(sit,area,version):
    """
    Calculate sea ice volume
    """
    sitq = sit / 1000.    
    
    volume = (sitq * area) / 1000.
    
    if version == 'PIOMAS':
        voln = np.squeeze(np.apply_over_axes(np.nansum,volume[:,:,:,:],(2,3)))
    else:
        voln = np.squeeze(np.apply_over_axes(np.nansum,volume[:,:,:,:],(3,4)))

    print 'Completed: Calculated volume!'    
    return voln


### Call functions   
sith,lat,lon = lens.readLENSEnsemble(directorydatal,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
sitp = readPIOMAS(directorydatap,0.15)
lons,lats = np.meshgrid(lon,lat)

area,lonsize,latsize = areaGrid(lat,lon)    
volh = volumeYear(sith,area,'historical')
volf = volumeYear(sitf,area,'rcp85')
volp = volumeYear(sitp,area,'PIOMAS')

#### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))
        else:
            spine.set_color('none')  
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([]) 

###########################################################################
###########################################################################
###########################################################################
###########################################################################
#volaveh = np.reshape(volh,(sith.shape[0],sith.shape[1]*sitf.shape[2]))
#volavef = np.reshape(volf,(sitf.shape[0],sitf.shape[1]*sitf.shape[2]))
#volavep = np.ravel(volp)  
#
#emptyf = np.empty((sitf.shape[0],1032))
#emptyf[:] = np.nan
#volavefq = np.append(emptyf,volavef,axis=1)
#
#emptyp = np.array([np.nan] * ((1979-1920)*12))
#newvolp = np.append(emptyp,volavep,axis=0)
#      
#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_linewidth(2)
#ax.spines['left'].set_linewidth(2)
#ax.tick_params('both',length=4.5,width=2,which='major')  
#
#for i in xrange(volaveh.shape[0]):
#    plt.plot(volaveh[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#for i in xrange(volavefq.shape[0]):
#    plt.plot(volavefq[i,:],color='forestgreen',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(newvolp,color='darkorange',alpha=1,linewidth=1,linestyle='-',
#         zorder=3,label=r'PIOMAS')
#         
#plt.axvline(1020,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,1933,240),np.arange(1920,2101,20))
#plt.yticks(np.arange(0,36,5),map(str,np.arange(0,36,5)))
#plt.xlim([0,1933])
#
#plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_vol_all.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
#volavehyr = np.nanmean(volh,axis=2)
#volavefyrq = np.nanmean(volf,axis=2)
#volavepyrq = np.nanmean(volp,axis=1)
#
#emptyf = np.empty((sitf.shape[0],sith.shape[1]))
#emptyf[:] = np.nan
#volavefyr = np.append(emptyf,volavefyrq,axis=1)
#
#emptyp = np.array([np.nan] * ((1979-1920)))
#volavepyr = np.append(emptyp,volavepyrq,axis=0)
#
#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_linewidth(2)
#ax.spines['left'].set_linewidth(2)
#ax.tick_params('both',length=4.5,width=2,which='major')  
#
#for i in xrange(volavehyr.shape[0]):
#    plt.plot(volavehyr[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#for i in xrange(volavefyr.shape[0]):
#    plt.plot(volavefyr[i,:],color='forestgreen',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(np.nanmean(volavehyr,axis=0),color='steelblue',linewidth=1.5,linestyle='-',
#         zorder=3,label=r'Mean Historical')    
#plt.plot(np.nanmean(volavefyr,axis=0),color='darkgreen',linewidth=1.5,linestyle='-',
#         zorder=3,label=r'Mean Future')   
#plt.plot(volavepyr,color='darkorange',alpha=1,linewidth=2,linestyle='-',
#         zorder=4,label=r'PIOMAS')
#         
#plt.axvline(85,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.yticks(np.arange(0,36,5),map(str,np.arange(0,36,5)))
#plt.xlim([0,180])
#
#plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_yr_vol_all.png',dpi=300)

############################################################################
############################################################################
############################################################################
#volcycle = np.nanmean(volh,axis=1)
#volcyclep = np.nanmean(volp,axis=0)
#
#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['bottom'].set_linewidth(2)
#ax.spines['left'].set_linewidth(2)
#ax.tick_params('both',length=4.5,width=2,which='major')  
##plt.grid(color='k',zorder=1,alpha=0.4)
#
#for i in xrange(volcycle.shape[0]):
#    plt.plot(volcycle[i],color='darkslateblue',alpha=0.3,linewidth=0.7,
#             zorder=2)
#plt.plot(np.nanmean(volcycle,axis=0),color='steelblue',linewidth=2,
#         marker='o',markeredgecolor='steelblue',zorder=3,
#         label=r'Mean LENS')
#         
#plt.plot(volcyclep,color='darkorange',linewidth=2,
#         marker='^',markeredgecolor='darkorange',
#         zorder=4,label=r'PIOMAS')
#
#plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')
#plt.xticks(np.arange(0,12,1),months) 
#plt.yticks(np.arange(0,36,5),map(str,np.arange(0,36,5)))
#plt.xlim([0,11]) 
#plt.ylim([0,35])
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'seasonalcycle_vol_lens.png',dpi=300)

############################################################################
############################################################################
############################################################################

volavehyrs = volh[:,:,8]
volavefyrqs = volf[:,:,8]
volavepyrqs = volp[:,8]

emptyf = np.empty((sitf.shape[0],sith.shape[1]))
emptyf[:] = np.nan
volavefyrs = np.append(emptyf,volavefyrqs,axis=1)

emptyp = np.array([np.nan] * ((1979-1920)))
volavepyrs = np.append(emptyp,volavepyrqs,axis=0)

volavehyrm = volh[:,:,2]
volavefyrqm = volf[:,:,2]
volavepyrqm = volp[:,2]

emptyf = np.empty((sitf.shape[0],sith.shape[1]))
emptyf[:] = np.nan
volavefyrm = np.append(emptyf,volavefyrqm,axis=1)

emptyp = np.array([np.nan] * ((1979-1920)))
volavepyrm = np.append(emptyp,volavepyrqm,axis=0)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(volavehyrs.shape[0]):
    plt.plot(volavehyrs[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
for i in xrange(volavefyrs.shape[0]):
    plt.plot(volavefyrs[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
for i in xrange(volavehyrs.shape[0]):
    plt.plot(volavehyrm[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
for i in xrange(volavefyrs.shape[0]):
    plt.plot(volavefyrm[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
plt.plot(np.nanmean(volavehyrs,axis=0),color='indianred',linewidth=1.5,linestyle='-',
         zorder=3,label=r'LENS September')  
plt.plot(np.nanmean(volavefyrs,axis=0),color='indianred',linewidth=1.5,linestyle='-',
         zorder=3)   
plt.plot(volavepyrs,color='darkorchid',alpha=1,linewidth=2,linestyle='-',
         zorder=4)
         
plt.plot(np.nanmean(volavehyrm,axis=0),color='teal',linewidth=1.5,linestyle='-',
         zorder=3,label=r'LENS March')    
plt.plot(np.nanmean(volavefyrm,axis=0),color='teal',linewidth=1.5,linestyle='-',
         zorder=3)   
plt.plot(volavepyrm,color='darkorchid',alpha=1,linewidth=2,linestyle='-',
         zorder=4,label=r'PIOMAS')
         
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
plt.yticks(np.arange(0,36,5),map(str,np.arange(0,36,5)))
plt.xlim([0,160])
plt.ylim([0,35])

plt.ylabel(r'\textbf{Sea Ice Volume ($\times$1000 km$^{3}$)}')

plt.legend(shadow=False,fontsize=9,loc='upper right',
           fancybox=True,frameon=False)

plt.savefig(directoryfigure + 'lens_sepmar_vol_all.png',dpi=300)

############################################################################
############################################################################
############################################################################
### Plot first ice free volume

#volavefyrqs = volf[:,:,8]
#
#zeros = np.empty((volavefyrqs.shape[0]))
#for i in xrange(volavefyrqs.shape[0]):
#    val = np.where(volavefyrqs[i,:] <= 1)[0]
#    zeros[i] = val[0] + 2006
#    
#fig = plt.figure()
#ax = plt.subplot(111)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#n, bins, patches = ax.hist(zeros,bins=range(2038,2056),align='left')
#
#for i in range(len(patches)):
#    patches[i].set_facecolor('cornflowerblue')
#    patches[i].set_edgecolor('white')
#    patches[i].set_linewidth(0.9)
#
#plt.xticks(np.arange(2038,2055,2),np.arange(2038,2055,2))
#plt.yticks(np.arange(0,7,1),map(str,np.arange(0,7,1)))
#plt.xlim([2038,2054])
#plt.ylim([0,6])
#
#plt.ylabel(r'\textbf{Number of Ensembles}')
#
#plt.savefig(directoryfigure + 'lens_vol_icefreetime_all.png',dpi=300)