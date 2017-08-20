"""
Script plots entire LENS time series for sea ice thickness during 
the months of March and September
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 5 June 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_LENS as lens
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
print '\n' '----LENS Mean Sea Ice Thickness - %s----' % titletime 

### Alott time series
yearmin = 1920
yearmax = 2080
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

### Call functions   
sith,lat,lon = lens.readLENSEnsemble(directorydatal,0.15,'historical')
sitf,lat,lon = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
sitp = readPIOMAS(directorydatap,0.15)
lons2,lats2 = np.meshgrid(lon,lat)

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
    
sitaveh = weightThick(sith,lats2,'lens')
sitavef = weightThick(sitf,lats2,'lens')
sitavep = weightThick(sitp,lats2,'piomas')

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

############################################################################
############################################################################
############################################################################

sitavehyrs = sitaveh[:,:,8]
sitavefyrqs = sitavef[:,:,8]
sitavepyrqs = sitavep[:,8]

emptyf = np.empty((sitf.shape[0],sith.shape[1]))
emptyf[:] = np.nan
sitavefyrs = np.append(emptyf,sitavefyrqs,axis=1)

emptyp = np.array([np.nan] * ((1979-1920)))
sitavepyrs = np.append(emptyp,sitavepyrqs,axis=0)

sitavehyrm = sitaveh[:,:,2]
sitavefyrqm = sitavef[:,:,2]
sitavepyrqm = sitavep[:,2]

emptyf = np.empty((sitf.shape[0],sith.shape[1]))
emptyf[:] = np.nan
sitavefyrm = np.append(emptyf,sitavefyrqm,axis=1)

emptyp = np.array([np.nan] * ((1979-1920)))
sitavepyrm = np.append(emptyp,sitavepyrqm,axis=0)

fig = plt.figure()
ax = plt.subplot(211)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitavehyrs.shape[0]):
    plt.plot(sitavehyrm[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
for i in xrange(sitavefyrs.shape[0]):
    plt.plot(sitavefyrm[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
        
         
plt.plot(np.nanmean(sitavehyrm,axis=0),color='teal',linewidth=1.5,linestyle='-',
         zorder=3)    
plt.plot(np.nanmean(sitavefyrm,axis=0),color='teal',linewidth=1.5,linestyle='-',
         zorder=3)   
plt.plot(sitavepyrm,color='darkorchid',alpha=1,linewidth=2,linestyle='-',
         zorder=4,label=r'PIOMAS')
         
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
plt.yticks(np.arange(0,5,1),map(str,np.arange(0,5,1))) 
plt.xlim([0,160])
plt.ylim([0,4])

plt.text(0,0,r'\textbf{MARCH}',fontsize=20,color='darkgrey')

plt.legend(shadow=False,fontsize=9,loc='upper right',
           fancybox=True,frameon=False)\

###########################################################################
###########################################################################
###########################################################################
ax = plt.subplot(212)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitavehyrs.shape[0]):
    plt.plot(sitavehyrs[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)        
for i in xrange(sitavefyrs.shape[0]):
    plt.plot(sitavefyrs[i,:],color='dimgrey',alpha=0.55,linewidth=0.3,
         zorder=2)
         
plt.plot(np.nanmean(sitavehyrs,axis=0),color='indianred',linewidth=1.5,linestyle='-',
         zorder=3,label=r'LENS September')  
plt.plot(np.nanmean(sitavefyrs,axis=0),color='indianred',linewidth=1.5,linestyle='-',
         zorder=3)   
plt.plot(sitavepyrs,color='darkorchid',alpha=1,linewidth=2,linestyle='-',
         zorder=4)

plt.axvline(85,linestyle='--',linewidth=2,color='k')
         
plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
plt.yticks(np.arange(0,5,1),map(str,np.arange(0,5,1))) 
plt.xlim([0,160])
plt.ylim([0,4])

plt.text(0,0,r'\textbf{SEPTEMBER}',fontsize=20,color='darkgrey')
plt.text(-18,6.6,r'\textbf{Sea Ice Thickness (m)}',fontsize=11,
                            rotation=90)
         
plt.subplots_adjust(hspace=0.3)

plt.savefig(directoryfigure + 'lens_sepmar_sit_all.png',dpi=300)