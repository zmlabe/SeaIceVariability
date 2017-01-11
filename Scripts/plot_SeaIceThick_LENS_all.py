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
print '\n' '----LENS Historical Mean Sea Ice Thickness - %s----' % titletime 

### Alott time series
yearmin = 1920
yearmax = 2100
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

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
#sitp = readPIOMAS(directorydatap,0.15)
#lons,lats = np.meshgrid(lons,lats)
#  
#sitaveh = weightThick(sith,lats,'lens')
#sitavef = weightThick(sitf,lats,'lens')
#sitavep = weightThick(sitp,lats,'piomas')

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
#sitaveh = np.reshape(sitaveh,(39,86*12))
#sitavef = np.reshape(sitavef,(39,75*12))
#sitavep = np.ravel(sitavep)  
#
#emptyf = np.empty((39,1032))
#emptyf[:] = np.nan
#sitavefq = np.append(emptyf,sitavef,axis=1)
#
#emptyp = np.array([np.nan] * ((1979-1920)*12))
#newsitp = np.append(emptyp,sitavep,axis=0)
#      
#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#plt.setp(ax.get_xticklabels(), visible=False)
#ax.xaxis.set_tick_params(size=0)
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#
#for i in xrange(sitaveh.shape[0]):
#    plt.plot(sitaveh[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#for i in xrange(sitavefq.shape[0]):
#    plt.plot(sitavefq[i,:],color='forestgreen',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(newsitp,color='darkorange',alpha=1,linewidth=1,linestyle='-',
#         zorder=3,label=r'PIOMAS')
#         
#plt.axvline(1020,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,1933,240),np.arange(1920,2101,20))
#plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
#plt.xlim([0,1933])
#
#plt.ylabel('Sea Ice Thickness (m)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_all.png',dpi=300)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
sitavehyr = np.nanmean(sitaveh,axis=2)
sitavefyrq = np.nanmean(sitavef,axis=2)
sitavepyrq = np.nanmean(sitavep,axis=1)

emptyf = np.empty((39,86))
emptyf[:] = np.nan
sitavefyr = np.append(emptyf,sitavefyrq,axis=1)

emptyp = np.array([np.nan] * ((1979-1920)))
sitavepyr = np.append(emptyp,sitavepyrq,axis=0)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitavehyr.shape[0]):
    plt.plot(sitavehyr[i,:],color='cornflowerblue',alpha=0.3,linewidth=0.7,
         zorder=2)
         
for i in xrange(sitavefyr.shape[0]):
    plt.plot(sitavefyr[i,:],color='forestgreen',alpha=0.3,linewidth=0.7,
         zorder=2)
         
plt.plot(np.nanmean(sitavehyr,axis=0),color='steelblue',linewidth=1.5,linestyle='-',
         zorder=3,label=r'Mean Historical')    
plt.plot(np.nanmean(sitavefyr,axis=0),color='darkgreen',linewidth=1.5,linestyle='-',
         zorder=3,label=r'Mean RCP85')   
plt.plot(sitavepyr,color='darkorange',alpha=1,linewidth=2,linestyle='-',
         zorder=4,label=r'PIOMAS')
         
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
plt.xlim([0,160])

plt.ylabel(r'\textbf{Sea Ice Thickness (m)}')

plt.legend(shadow=False,fontsize=9,loc='upper right',
           fancybox=True,frameon=False)

plt.savefig(directoryfigure + 'lens_yr_all.png',dpi=300)

############################################################################
############################################################################
############################################################################