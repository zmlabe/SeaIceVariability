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
from mpl_toolkits.basemap import Basemap
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
sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
    
sitp = readPIOMAS(directorydatap,0.15)
lons2,lats2 = np.meshgrid(lons,lats)

### Slice regions
def RegionalSlice(sith,sitf,sitp,region,lats,lons):
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

    latq = np.where((lats >= latmin) & (lats <= latmax))[0]
    latsn = lats[latq]
    lonq = np.where((lons >= lonmin) & (lons <= lonmax))[0]
    lonsn = lons[lonq]
    
    lons2n,lats2n = np.meshgrid(lonsn,latsn)

    # Slice lats
    sitmh = sith[:,:,:,latq,:]
    sitmf = sitf[:,:,:,latq,:]
    sitmp = sitp[:,:,latq,:]
    
    # Slice lons
    sitmmh = sitmh[:,:,:,:,lonq]
    sitmmf = sitmf[:,:,:,:,lonq]
    sitmmp = sitmp[:,:,:,lonq]
    
    sitall = np.append(sitmmh,sitmmf,axis=1)
    
    return sitall,sitmmp,lats2n,lons2n
    
### Regions sliced   
sitallg,sitppg,lats2ng,lons2ng = RegionalSlice(sith,sitf,sitp,'Greenland',lats,lons)
sitallb,sitppb,lats2nb,lons2nb = RegionalSlice(sith,sitf,sitp,'BeaufortSea',lats,lons)
sitalle,sitppe,lats2ne,lons2ne = RegionalSlice(sith,sitf,sitp,'EastSiberianSea',lats,lons)
sitalll,sitppl,lats2nl,lons2nl = RegionalSlice(sith,sitf,sitp,'LaptevSea',lats,lons)
sitallkb,sitppkb,lats2nkb,lons2nkb = RegionalSlice(sith,sitf,sitp,'KB',lats,lons)
sitallcab,sitppcab,lats2ncab,lons2ncab = RegionalSlice(sith,sitf,sitp,'CAB',lats,lons)

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

### Calculate time series per region  
sityrg = weightThick(sitallg,lats2ng,'lens')    
sityrb = weightThick(sitallb,lats2nb,'lens') 
sityre = weightThick(sitalle,lats2ne,'lens') 
sityrl = weightThick(sitalll,lats2nl,'lens') 
sityrkb = weightThick(sitallkb,lats2nkb,'lens') 
sityrcab = weightThick(sitallcab,lats2ncab,'lens')   

sitppyrg = weightThick(sitppg,lats2ng,'piomas')    
sitppyrb = weightThick(sitppb,lats2nb,'piomas') 
sitppyre = weightThick(sitppe,lats2ne,'piomas') 
sitppyrl = weightThick(sitppl,lats2nl,'piomas') 
sitppyrkb = weightThick(sitppkb,lats2nkb,'piomas') 
sitppyrcab = weightThick(sitppcab,lats2ncab,'piomas')     

### Select month
monthq = 8
if monthq < 12:
    sityrg = sityrg[:,:,monthq]
    sityrb = sityrb[:,:,monthq]
    sityre = sityre[:,:,monthq]
    sityrl = sityrl[:,:,monthq]
    sityrkb = sityrkb[:,:,monthq]
    sityrcab = sityrcab[:,:,monthq]
    
    sitppyrg = sitppyrg[:,monthq]
    sitppyrb = sitppyrb[:,monthq]
    sitppyre = sitppyre[:,monthq]
    sitppyrl = sitppyrl[:,monthq]
    sitppyrkb = sitppyrkb[:,monthq]
    sitppyrcab = sitppyrcab[:,monthq]
    
yearq = 'on'
if yearq == 'on':
    yearq = np.where((years >= 1979) & (years <= 2015))[0]
    
    sityrg = sityrg[:,yearq]
    sityrb = sityrb[:,yearq]
    sityre = sityre[:,yearq]
    sityrl = sityrl[:,yearq]
    sityrkb = sityrkb[:,yearq]
    sityrcab = sityrcab[:,yearq]
    
### Calculate mean per region
sitmeang = np.nanmean(sityrg,axis=0)  
sitmeanb = np.nanmean(sityrb,axis=0)  
sitmeane = np.nanmean(sityre,axis=0)  
sitmeanl = np.nanmean(sityrl,axis=0)   
sitmeankb = np.nanmean(sityrkb,axis=0)  
sitmeancab = np.nanmean(sityrcab,axis=0)   
        
#### Plot Figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
### Adjust axes in time series plots 
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
        
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
         
plt.plot(sitmeang,color='steelblue',linewidth=2,linestyle='-',
         zorder=3,label=r'Greenland Sea')    
plt.plot(sitmeanb,color='darkgreen',linewidth=2,linestyle='-',
         zorder=4,label=r'Beaufort Sea')   
plt.plot(sitmeane,color='darkorange',linewidth=2,linestyle='-',
         zorder=5,label=r'East Siberian Sea')
plt.plot(sitmeanl,color='darkblue',linewidth=2,linestyle='-',
         zorder=6,label=r'Laptev Sea')    
plt.plot(sitmeankb,color='m',linewidth=2,linestyle='-',
         zorder=7,label=r'Kara-Barents Seas')   
plt.plot(sitmeancab,color='yellowgreen',linewidth=2,linestyle='-',
         zorder=8,label=r'Central Arctic Basin')
         
plt.plot(sitppyrg,color='steelblue',linewidth=1,linestyle='--',
         zorder=3)    
plt.plot(sitppyrb,color='darkgreen',linewidth=1,linestyle='--',
         zorder=4)   
plt.plot(sitppyre,color='darkorange',linewidth=1,linestyle='--',
         zorder=5)
plt.plot(sitppyrl,color='darkblue',linewidth=1,linestyle='--',
         zorder=6)    
plt.plot(sitppyrkb,color='m',linewidth=1,linestyle='--',
         zorder=7)   
plt.plot(sitppyrcab,color='yellowgreen',linewidth=1,linestyle='--',
         zorder=8)

xlabels = map(str,np.arange(1979,2017,5))
plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.yticks(np.arange(0,6,1),map(str,np.arange(0,6,1))) 
plt.ylim([0,5])

plt.ylabel(r'\textbf{Sea Ice Thickness (m)}')

plt.legend(shadow=False,fontsize=7,loc='upper right',
           fancybox=True,frameon=False,ncol=1)
           
plt.savefig(directoryfigure+'lens_yr_regional_sept_sat.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.plot([0]*len(sitmeang),color='k',linewidth=3,linestyle='--')
         
plt.plot(sitmeang - sitppyrg,color='steelblue',linewidth=2,linestyle='-',
         zorder=3,label=r'Greenland Sea')    
plt.plot(sitmeanb - sitppyrb,color='darkgreen',linewidth=2,linestyle='-',
         zorder=4,label=r'Beaufort Sea')   
plt.plot(sitmeane - sitppyre,color='darkorange',linewidth=2,linestyle='-',
         zorder=5,label=r'East Siberian Sea')
plt.plot(sitmeanl - sitppyrl,color='darkblue',linewidth=2,linestyle='-',
         zorder=6,label=r'Laptev Sea')    
plt.plot(sitmeankb - sitppyrkb,color='m',linewidth=2,linestyle='-',
         zorder=7,label=r'Kara-Barents Seas')   
plt.plot(sitmeancab - sitppyrcab,color='yellowgreen',linewidth=2,linestyle='-',
         zorder=8,label=r'Central Arctic Basin')
         
xlabels = map(str,np.arange(1979,2017,5))
plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-3,3])

plt.ylabel(r'\textbf{Sea Ice Thickness (m)}')

plt.legend(shadow=False,fontsize=7,loc='bottom right',
           fancybox=True,frameon=False,ncol=1)
                      
plt.savefig(directoryfigure+'lens_yr_regional_sep_sat_difference.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################

fig = plt.figure(figsize=(5,9))

ax = plt.subplot(612)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=0,width=0,which='major',color='none')
ax.xaxis.set_visible(False)

plt.plot(sityrg-sitppyrg,color='dimgrey',alpha=0.55,linewidth=0.2)
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeang - sitppyrg,color='steelblue',linewidth=2,linestyle='-',
         zorder=3,label=r'Greenland Sea',marker='o',markersize=4,
         markeredgecolor='steelblue')   
         
#xlabels = map(str,np.arange(1979,2017,5))
#plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{GREENLAND}',fontsize=11,color='darkgrey')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])
         
###########################################################################         
         
ax = plt.subplot(613)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=0,width=0,which='major',color='none')
ax.xaxis.set_visible(False)

plt.plot(sityrb-sitppyrb,color='dimgrey',alpha=0.55,linewidth=0.2)
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeanb - sitppyrb,color='darkgreen',linewidth=2,linestyle='-',
         zorder=4,label=r'Beaufort Sea',marker='o',markersize=4,
         markeredgecolor='darkgreen')

#xlabels = map(str,np.arange(1979,2017,5))
#plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{BEAUFORT-CHUKCHI SEAS}',fontsize=11,color='darkgrey')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])
         
###########################################################################          
         
ax = plt.subplot(614)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=0,width=0,which='major',color='none')
ax.xaxis.set_visible(False)

plt.plot(sityre-sitppyre,color='dimgrey',alpha=0.55,linewidth=0.2)         
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeane - sitppyre,color='darkorange',linewidth=2,linestyle='-',
         zorder=5,label=r'East Siberian Sea',marker='o',markersize=4,
         markeredgecolor='darkorange')

#xlabels = map(str,np.arange(1979,2017,5))
#plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{EAST SIBERIAN SEA}',fontsize=11,color='darkgrey')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])

########################################################################### 

ax = plt.subplot(615)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=0,width=0,which='major',color='none')
ax.xaxis.set_visible(False)

plt.plot(sityrl-sitppyrl,color='dimgrey',alpha=0.55,linewidth=0.2)         
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeanl - sitppyrl,color='darkblue',linewidth=2,linestyle='-',
         zorder=6,label=r'Laptev Sea',marker='o',markersize=4,
         markeredgecolor='darkblue')   

#xlabels = map(str,np.arange(1979,2017,5))
#plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{LAPTEV SEA}',fontsize=11,color='darkgrey')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])

########################################################################### 

ax = plt.subplot(616)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=4,width=1.5,which='major',color='darkgrey')
 
plt.plot(sityrkb-sitppyrkb,color='dimgrey',alpha=0.55,linewidth=0.2)        
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeankb - sitppyrkb,color='m',linewidth=2,linestyle='-',
         zorder=7,label=r'Kara-Barents Seas',marker='o',markersize=4,
         markeredgecolor='m')   

xlabels = map(str,np.arange(1979,2017,5))
plt.xticks(np.arange(0,39,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{BARENTS-KARA SEAS}',fontsize=11,color='darkgrey')
plt.text(-5.35,15.25,r'\textbf{Difference (m)}',fontsize=15,
                          color='k',rotation=90)

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])

########################################################################### 

ax = plt.subplot(611)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('none')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params('x',length=0,width=0,which='major',color='none')
ax.xaxis.set_visible(False)

plt.plot(sityrcab-sitppyrcab,color='dimgrey',alpha=0.55,linewidth=0.2)
plt.plot([0]*len(sitmeang),color='k',linewidth=1.5,linestyle='--')
plt.plot(sitmeancab - sitppyrcab,color='r',linewidth=2,linestyle='-',
         zorder=8,label=r'Central Arctic Basin',marker='o',markersize=4,
         markeredgecolor='r')
         
#xlabels = map(str,np.arange(1979,2006,5))
#plt.xticks(np.arange(0,28,5),xlabels)
plt.xlim([0,36])

plt.text(0,2.8,r'\textbf{CENTRAL ARCTIC BASIN}',fontsize=11,color='darkgrey')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1))) 
plt.ylim([-1,3])

plt.savefig(directoryfigure + 'SIT_regions_diff_sat_sep.png',dpi=300)