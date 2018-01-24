"""
Scripts calculates SIT trends from LENS
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 23 February 2017
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
import scipy.stats as sts

### Define directories
directorydatal = '/surtsey/ypeings/'
directorydatap = '/surtsey/zlabe/seaice_obs/PIOMAS/Thickness/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'
directorydata2 = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----LENS Historical Mean Sea Ice Thickness - %s----' % titletime 

## Alott time series
yearmin = 1920
yearmax = 2080
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))
    
### Call functions   
sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
lons2,lats2 = np.meshgrid(lons,lats)

### Slice regions
def RegionalSlice(sith,sitf,region,lats,lons):
    """
    Mask out selected marginal seas and their time series
    """

    if region == 'BeaufortSea':
        latmin = 68
        latmax = 85
        lonmin = 185
        lonmax = 235
    elif region == 'Greenland':
        latmin = 76
        latmax = 85
        lonmin = 235
        lonmax = 370
        lonmin2 = 0
        lonmax2 = 10
    elif region == 'EastSiberianSea':
        latmin = 68
        latmax = 85
        lonmin = 146
        lonmax = 185
    elif region == 'LaptevSea':
        latmin = 70
        latmax = 85
        lonmin = 100
        lonmax = 146
    elif region == 'KaraSea':
        latmin = 68
        latmax = 85
        lonmin = 50
        lonmax = 100
    elif region == 'BarentsSea':
        latmin = 68
        latmax = 85
        lonmin = 10
        lonmax = 50
    elif region == 'CAB':
        latmin = 85
        latmax = 90
        lonmin = 0
        lonmax = 360
    elif region == 'KB':
        latmin = 68
        latmax = 85
        lonmin = 10
        lonmax = 100
    
    if region == 'Greenland':
        latq = np.where((lats >= latmin) & (lats <= latmax))[0]
        latsn = lats[latq]
        lonq1 = np.where((lons >= lonmin) & (lons <= lonmax))[0]
        lonq2 = np.where((lons >= lonmin2) & (lons <= lonmax2))[0]
        lonq = np.append(lonq2,lonq1,axis=0)
        lonsn = lons[lonq]
        lons2n,lats2n = np.meshgrid(lonsn,latsn)
        
    else:
        latq = np.where((lats >= latmin) & (lats <= latmax))[0]
        latsn = lats[latq]
        lonq = np.where((lons >= lonmin) & (lons <= lonmax))[0]
        lonsn = lons[lonq]
        lons2n,lats2n = np.meshgrid(lonsn,latsn)

    # Slice lats
    sitmh = sith[:,:,:,latq,:]
    sitmf = sitf[:,:,:,latq,:]
    
    # Slice lons
    sitmmh = sitmh[:,:,:,:,lonq]
    sitmmf = sitmf[:,:,:,:,lonq]
    
    sitall = np.append(sitmmh,sitmmf,axis=1)
    
    return sitall,lats2n,lons2n
    
#### Regions sliced   
sitallg,lats2ng,lons2ng = RegionalSlice(sith,sitf,'Greenland',lats,lons)
sitallb,lats2nb,lons2nb = RegionalSlice(sith,sitf,'BeaufortSea',lats,lons)
sitalle,lats2ne,lons2ne = RegionalSlice(sith,sitf,'EastSiberianSea',lats,lons)
sitalll,lats2nl,lons2nl = RegionalSlice(sith,sitf,'LaptevSea',lats,lons)
sitallkb,lats2nkb,lons2nkb = RegionalSlice(sith,sitf,'KB',lats,lons)
sitallcab,lats2ncab,lons2ncab = RegionalSlice(sith,sitf,'CAB',lats,lons)

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
        varq = var[:,:]
        mask = np.isfinite(varq) & np.isfinite(lats)
        varmask = varq[mask]
        areamask = np.cos(np.deg2rad(lats[mask]))
        sityr = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr

### Calculate time series per region  
sityrg = weightThick(sitallg,lats2ng,'lens')    
sityrb = weightThick(sitallb,lats2nb,'lens') 
sityre = weightThick(sitalle,lats2ne,'lens') 
sityrl = weightThick(sitalll,lats2nl,'lens') 
sityrkb = weightThick(sitallkb,lats2nkb,'lens') 
sityrcab = weightThick(sitallcab,lats2ncab,'lens')     

##### Select month (mar/sep)
monthqs = 8
sityrgs = sityrg[:,:,monthqs]
sityrbs = sityrb[:,:,monthqs]
sityres = sityre[:,:,monthqs]
sityrls = sityrl[:,:,monthqs]
sityrkbs = sityrkb[:,:,monthqs]
sityrcabs = sityrcab[:,:,monthqs]

#sitmeang = np.nanmean(sityrg,axis=0)  
#sitmeanb = np.nanmean(sityrb,axis=0)  
#sitmeane = np.nanmean(sityre,axis=0)  
#sitmeanl = np.nanmean(sityrl,axis=0)   
#sitmeankb = np.nanmean(sityrkb,axis=0)  
#sitmeancab = np.nanmean(sityrcab,axis=0)   

#sityrgs = np.nanmean(sityrg[:,:,:],axis=2)
#sityrbs = np.nanmean(sityrb[:,:,:],axis=2)
#sityres = np.nanmean(sityre[:,:,:],axis=2)
#sityrls = np.nanmean(sityrl[:,:,:],axis=2)
#sityrkbs = np.nanmean(sityrkb[:,:,:],axis=2)
#sityrcabs = np.nanmean(sityrcab[:,:,:],axis=2)

### Calculate first timing when threshold is reached 
thresh = 0.5
zerosgs = np.empty((sityrgs.shape[0]))
zerosbs = np.empty((sityrbs.shape[0]))
zeroses = np.empty((sityres.shape[0]))
zerosls = np.empty((sityrls.shape[0]))
zeroskbs = np.empty((sityrkbs.shape[0]))
zeroscabs = np.empty((sityrcabs.shape[0]))
for i in xrange(sityrgs.shape[0]):
    valgs = np.where(sityrgs[i,:] <= thresh)[0]
    zerosgs[i] = valgs[0]    
    vales = np.where(sityres[i,:] <= thresh)[0]
    zeroses[i] = vales[0] 
    valbs = np.where(sityrbs[i,:] <= thresh)[0]
    zerosbs[i] = valbs[0] 
    valls = np.where(sityrls[i,:] <= thresh)[0]
    zerosls[i] = valls[0] 
    valkbs = np.where(sityrkbs[i,:] <= thresh)[0]
    zeroskbs[i] = valkbs[0] 
    valcabs = np.where(sityrcabs[i,:] <= thresh)[0]
    zeroscabs[i] = valcabs[0] 
    

#valgs = np.where(sitmeang[:] <= thresh)[0]
#zerosgmean = years[valgs[0]]    
#vales = np.where(sitmeane[:] <= thresh)[0]
#zerosemean = years[vales[0]] 
#valbs = np.where(sitmeanb[:] <= thresh)[0]
#zerosbmean = years[valbs[0]] 
#valls = np.where(sitmeanl[:] <= thresh)[0]
#zeroslmean = years[valls[0]] 
#valkbs = np.where(sitmeankb[:] <= thresh)[0]
#zeroskbmean = years[valkbs[0]]
#valcabs = np.where(sitmeancab[:] <= thresh)[0]
#zeroscabmean = years[valcabs[0]] 
    
sep = [zeroscabs,zerosgs,zerosbs,zeroses,zerosls,zeroskbs]
yearsep = []
for i in range(len(sep)):
    yearsepq = years[sep[i].astype(int)]
    yearsep.append(yearsepq)
    
meansep = []
perc05sep = []
perc95sep = []
for i in range(len(sep)):
    meansepq = np.mean(yearsep[i])
    perc05sepq = np.percentile(yearsep[i],5)
    perc95sepq = np.percentile(yearsep[i],95)
    
    meansep.append(meansepq)
    perc05sep.append(perc05sepq)
    perc95sep.append(perc95sepq)
    
meansep = np.asarray(meansep)
perc05sep = np.asarray(perc05sep)
perc95sep = np.asarray(perc95sep)
        
#### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
### Adjust axes in time series plots 
def adjust_spines(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 0))
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('dimgrey')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4,width=2,which='major',color='dimgrey',
               labelbottom='off',bottom='off')

#ax.yaxis.grid(zorder=0,color='darkgrey',alpha=0.35,linestyle='-')

xlabels = [r'\textbf{CENTRAL ARCTIC',r'\textbf{GREENLAND}',
           r'\textbf{BEAUFORT-CHUKCHI}',r'\textbf{EAST-SIBERIAN}',
           r'\textbf{LAPTEV}',r'\textbf{BARENTS-KARA}']
ccc= ['r','steelblue','darkgreen','darkorange','darkblue','m']
for i in range(len(meansep)):
    plt.scatter(i,meansep[i],s=100,c=ccc[i],edgecolor=ccc[i],zorder=5)
    plt.errorbar(i,meansep[i],
                 yerr=np.array([[meansep[i]-perc05sep[i],perc95sep[i]-meansep[i]]]).T,
                 color=ccc[i],linewidth=1.5,capthick=3,capsize=10)
    print([meansep[i]-perc05sep[i],perc95sep[i]-meansep[i]])
    print([perc05sep[i],perc95sep[i]])
    
    plt.text(i,perc95sep[i]+3,r'\textbf{%s}' % xlabels[i],
             color='dimgrey',fontsize=9,ha='center',
             va='center')

plt.xticks(np.arange(-1,7,1),xlabels,rotation=15,fontsize=7)
plt.yticks(np.arange(1990,2080,10),list(map(str,np.arange(1990,2080,10))))

plt.ylim([1995,2070])
plt.xlim([-1,6])

plt.savefig(directoryfigure + 'timing_SIT_LENS_regional.png',dpi=300)
         
