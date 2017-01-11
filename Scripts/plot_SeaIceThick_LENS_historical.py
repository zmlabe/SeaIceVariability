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
yearmax = 2005
years = np.arange(yearmin,yearmax+1,1)
yearf = np.arange(2006,2080+1,1)

yearp = np.arange(1979,2015+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
#ensemble = np.array(['02','03','04','05','06','07','08','09'])
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

#### Call functions
#sit,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
#lons,lats = np.meshgrid(lons,lats)
#
#sitp = readPIOMAS(directorydatap,0.15)

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

#### Call functions     
#sitave = weightThick(sit,lats,'lens')
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
#sitave = np.reshape(sitave,(39,86*12))
#sitavep = np.ravel(sitavep)  
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
#ax.spines['bottom'].set_linewidth(2)
#ax.spines['left'].set_linewidth(2)
#ax.tick_params('both',length=4.5,width=2,which='major')  
#
#for i in xrange(sitave.shape[0]):
#    plt.plot(sitave[i,:],color='darkslateblue',alpha=0.3,linewidth=0.7,
#         zorder=2)
#         
#plt.plot(newsitp,color='darkorange',alpha=1,linewidth=1,linestyle='-',
#         zorder=3,label=r'PIOMAS')
#
#plt.xticks(np.arange(0,1161,120),np.arange(1920,2011,10))
#plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
#plt.xlim([0,1080])
#
#plt.ylabel('Sea Ice Thickness (m)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_historical.png',dpi=300)


###########################################################################
###########################################################################
###########################################################################
###########################################################################
#sitaveyr = np.nanmean(sitave,axis=2)
#sitavepyr = np.nanmean(sitavep,axis=1)
#
#fig = plt.figure()
#ax = plt.subplot(111)
#
#emptyp = np.array([np.nan] * ((1979-1920)))
#newsitp = np.append(emptyp,sitavepyr,axis=0)
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
#
#plt.plot(np.nanmean(sitaveyr,axis=0),color='steelblue',linewidth=1.5,linestyle='-',
#         zorder=3,label=r'Mean LENS')         
#plt.plot(newsitp,color='darkorange',alpha=1,linewidth=2,linestyle='-',
#         zorder=4,label=r'PIOMAS')
#
#plt.xticks(np.arange(0,96,10),np.arange(1920,2011,10))
#plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
#plt.xlim([0,95])
#
#plt.ylabel('Sea Ice Thickness (m)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper right',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_historical_yr.png',dpi=300)

############################################################################
############################################################################
############################################################################
yearlq = np.where((years >= 1979) & (years <= 2005))[0]
yearpq = np.where((yearp >= 1979) & (yearp <= 2005))[0]
sitcycle = np.nanmean(sitave[:,yearlq,:],axis=1)
sitcyclep = np.nanmean(sitavep[yearpq,:],axis=0)

yearlqf = np.where((yearf >= 2006) & (yearf <= 2015))[0]
yearpqf = np.where((yearp >= 2006) & (yearp <= 2015))[0]
sitcyclef = np.nanmean(sitavef[:,yearlqf,:],axis=1)
sitcyclefp = np.nanmean(sitavep[yearpqf,:],axis=0)

fig = plt.figure()
ax = plt.subplot(211)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
plt.setp(ax.get_xticklabels(), visible=False)
ax.xaxis.set_tick_params(size=0)
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitcycle.shape[0]):
    plt.plot(sitcycle[i],color='cornflowerblue',alpha=0.3,linewidth=0.7,
             zorder=2)
plt.plot(np.nanmean(sitcycle,axis=0),color='darkorchid',linewidth=2,
         marker='o',markeredgecolor='darkorchid',zorder=3,
         label=r'Mean LENS',markersize=4)
         
plt.plot(sitcyclep,color='seagreen',linewidth=2,
         marker='^',markeredgecolor='seagreen',
         zorder=4,label=r'PIOMAS',markersize=4)

plt.xticks(np.arange(0,12,1),months) 
plt.yticks(np.arange(0,5,1),map(str,np.arange(0,5,1))) 
plt.xlim([0,11]) 
plt.ylim([1,4])

plt.legend(shadow=False,fontsize=9,loc='upper right',
           fancybox=True,frameon=False)
           
plt.annotate(r'\textbf{HISTORICAL}', xy=(0, 0), xytext=(0.01,0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')           
           
ax = plt.subplot(212)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitcycle.shape[0]):
    plt.plot(sitcyclef[i],color='cornflowerblue',alpha=0.3,linewidth=0.7,
             zorder=2)
plt.plot(np.nanmean(sitcyclef,axis=0),color='darkorchid',linewidth=2,
         marker='o',markeredgecolor='darkorchid',zorder=3,
         label=r'Mean LENS',markersize=4)
         
plt.plot(sitcyclefp,color='seagreen',linewidth=2,
         marker='^',markeredgecolor='seagreen',
         zorder=4,label=r'PIOMAS',markersize=4)

plt.xticks(np.arange(0,12,1),months) 
plt.yticks(np.arange(0,5,1),map(str,np.arange(0,5,1))) 
plt.xlim([0,11]) 
plt.ylim([1,4])   

plt.annotate(r'\textbf{RCP8.5}', xy=(0, 0), xytext=(0.01,0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')  
            
plt.annotate(r'\textbf{Sea Ice Thickness (m)}', xy=(0, 0), xytext=(0.03,0.75),
            xycoords='figure fraction',fontsize=15,color='k',
            rotation=90)               

plt.savefig(directoryfigure + 'seasonalcycle_lens.png',dpi=300)