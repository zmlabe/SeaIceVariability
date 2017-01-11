"""
Scripts plots SST for CESM Large ensemble over historical and future
periods

Notes
-----
    Reference : Kay et al. [2014]
    Author : Zachary Labe
    Date   : 1 December 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_var_LENS as LV

### Define directories
directorydata = '/home/zlabe/Surtsey3/CESM_large_ensemble/' 
directoryfigure = '/home/zlabe/Desktop/'

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

sstk[np.where(sstk < 0)] = np.nan

## Convert to celsius 
sst = sstk - 273.15 

### Mask freezing
sst[np.where(sst <= -1.8)] = -1.8 # freezing point SST-ice

### Meshgrid of lats/lons
lons2,lats2 = np.meshgrid(lons,lats)

def weightAve(var,lats,types):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    """
    
    if types == 'lens':
        varyr = np.empty((var.shape[0],var.shape[1],var.shape[2]))
        for ens in xrange(var.shape[0]):
            for i in xrange(var.shape[1]):
                for j in xrange(var.shape[2]):
                    varq = var[ens,i,j,:,:]
                    mask = np.isfinite(varq) & np.isfinite(lats)
                    varmask = varq[mask]
                    areamask = np.cos(np.deg2rad(lats[mask]))
                    varyr[ens,i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
            
            print 'Completed: Weighting per ensemble #%s!' % ensembles[ens]
     
    print '\nCompleted: Yearly weighted SST average!' 
    return varyr
    
sstave = weightAve(sst,lats2,'lens')

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
#sstaven = np.nanmean(sstave,axis=2)
#sstaveq = np.nanmean(sstaven,axis=0)
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
#plt.plot(sstaven.transpose(),color='cornflowerblue',alpha=0.3,linewidth=0.7,
#     zorder=2)
#plt.plot(sstaveq,color='darkorchid',alpha=1,linewidth=2,
#     zorder=2)
#         
#plt.axvline(85,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.yticks(np.arange(-2,2.1,0.5),map(str,np.arange(-2,2.1,0.5)))
#plt.xlim([0,180])
#
#plt.ylabel(r'Sea Surface Temperatures ($^{\circ}$C)')
#
#plt.savefig(directoryfigure + 'lens_sstmean70N_all.png',dpi=300)
#
############################################################################
############################################################################
############################################################################
############################################################################
#sstsep = sstave[:,:,8]
#sstmar = sstave[:,:,2]
#
#sstsepmean = np.nanmean(sstave[:,:,8],axis=0)
#sstmarmean = np.nanmean(sstave[:,:,2],axis=0)
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
#plt.plot(sstsep.transpose(),color='darkgrey',alpha=0.55,linewidth=0.7,
#     zorder=2)
#plt.plot(sstsepmean,color='indianred',alpha=1,linewidth=2,
#     zorder=2,label=r'LENS September')
#     
#plt.plot(sstmar.transpose(),color='darkgrey',alpha=0.55,linewidth=0.7,
#     zorder=2)
#plt.plot(sstmarmean,color='steelblue',alpha=1,linewidth=2,
#     zorder=2,label=r'LENS March')
#     
#plt.axvline(85,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.yticks(np.arange(-2,8,1),map(str,np.arange(-2,8,1)))
#plt.xlim([0,180])
#
#plt.ylabel(r'Sea Surface Temperatures ($^{\circ}$C)')
#
#plt.legend(shadow=False,fontsize=11,loc='upper left',
#           fancybox=True,frameon=False)
#
#plt.savefig(directoryfigure + 'lens_sstmean70N_all_marsep.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
###########################################################################
sst_w = np.nanmean(sstave[:,:,0:3],axis=2)
sst_sp = np.nanmean(sstave[:,:,3:6],axis=2)
sst_su = np.nanmean(sstave[:,:,6:9],axis=2)
sst_f = np.nanmean(sstave[:,:,9:12],axis=2)

sstmean_w = np.nanmean(sst_w[:,:],axis=0)
sstmean_sp = np.nanmean(sst_sp[:,:],axis=0)
sstmean_su = np.nanmean(sst_su[:,:],axis=0)
sstmean_f = np.nanmean(sst_f[:,:],axis=0)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.plot(sst_w.transpose(),color='dimgrey',alpha=0.55,linewidth=0.3,
     zorder=2)
plt.plot(sstmean_w,color='darkslateblue',alpha=1,linewidth=2,
     zorder=4,label=r'JFM')
     
plt.plot(sst_sp.transpose(),color='dimgrey',alpha=0.55,linewidth=0.3,
     zorder=2)
plt.plot(sstmean_sp,color='teal',alpha=1,linewidth=2,
     zorder=5,label=r'AMJ')
     
plt.plot(sst_su.transpose(),color='dimgrey',alpha=0.55,linewidth=0.3,
     zorder=2)
plt.plot(sstmean_su,color='seagreen',alpha=1,linewidth=2,
     zorder=6,label=r'JAS')
     
plt.plot(sst_f.transpose(),color='dimgrey',alpha=0.55,linewidth=0.3,
     zorder=2)
plt.plot(sstmean_f,color='saddlebrown',alpha=1,linewidth=2,
     zorder=7,label=r'OND')
     
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
plt.yticks(np.arange(-2,8,1),map(str,np.arange(-2,8,1)))
plt.xlim([0,160])

plt.ylabel(r'\textbf{Sea Surface Temperatures ($^{\circ}$C)}')

plt.legend(shadow=False,fontsize=9,loc='upper left',
           fancybox=True,frameon=False)

plt.savefig(directoryfigure + 'lens_sstmean70N_all_seasons.png',dpi=300)
