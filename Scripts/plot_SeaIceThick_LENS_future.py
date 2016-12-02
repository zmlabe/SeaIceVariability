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
directoryfigure = '/home/zlabe/Desktop/LENS/'
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
yearmax = 2100
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + map(str,np.arange(10,16))
#ensemble = ['02','03','04','05','06','07','08','09'] + \
#        map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

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

#### Call functions
sit,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
lons,lats = np.meshgrid(lons,lats)

sitp = readPIOMAS(directorydatap,0.15)

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
sitave = weightThick(sit,lats,'lens')
sitavep = weightThick(sitp,lats,'piomas')

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
#sitave = np.reshape(sitave,(14,75*12))
#sitavep = np.ravel(sitavep)  
#
#timercp85 = len(np.where((years <= 2015))[0])*12
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
#for i in xrange(sitave.shape[0]):
#    plt.plot(sitave[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(sitavep[-(timercp85+1):],color='darkorange',alpha=1,linewidth=1,linestyle='-',
#         zorder=3,label=r'PIOMAS')
#
#plt.xticks(np.arange(0,1161,120),np.arange(2006,2101,10))
#plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
#plt.xlim([0,1080])
#
#plt.ylabel('Sea Ice Thickness (m)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_rcp85.png',dpi=300)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
#sitaveyr = np.nanmean(sitave,axis=2)
#sitavepyr = np.nanmean(sitavep,axis=1)
#
#timercp85 = len(np.where((years <= 2015))[0])
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
#for i in xrange(sitave.shape[0]):
#    plt.plot(sitaveyr[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(np.nanmean(sitaveyr,axis=0),color='steelblue',linewidth=1.5,linestyle='-',
#         zorder=3,label=r'Mean LENS')         
#plt.plot(sitavepyr[-(timercp85+1):],color='darkorange',alpha=1,linewidth=2,linestyle='-',
#         zorder=4,label=r'PIOMAS')
#
#plt.xticks(np.arange(0,96,10),np.arange(2006,2101,10))
#plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
#plt.xlim([0,95])
#
#plt.ylabel('Sea Ice Thickness (m)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_yr_rcp85.png',dpi=300)

############################################################################
############################################################################
############################################################################
sitcycle = np.nanmean(sitave,axis=1)
sitcyclep = np.nanmean(sitavep[-11:],axis=0)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.5,width=2,which='major')  
#plt.grid(color='k',zorder=1,alpha=0.4)

for i in xrange(sitcycle.shape[0]):
    plt.plot(sitcycle[i],color='darkslateblue',alpha=0.3,linewidth=0.7,
             zorder=2)
plt.plot(np.nanmean(sitcycle,axis=0),color='steelblue',linewidth=2,
         marker='o',markeredgecolor='steelblue',zorder=3,
         label=r'Mean LENS')
         
plt.plot(sitcyclep,color='darkorange',linewidth=2,
         marker='^',markeredgecolor='darkorange',
         zorder=4,label=r'PIOMAS')

plt.ylabel('Sea Ice Thickness (m)')
plt.xticks(np.arange(0,12,1),months) 
plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
plt.xlim([0,11]) 
plt.ylim([0,4])

plt.legend(shadow=False,fontsize=11,loc='upper right',
           fancybox=True,frameon=False)

plt.savefig(directoryfigure + 'seasonalcycle_lens_rcp85.png',dpi=300)