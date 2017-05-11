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
directorydatal = '/home/zlabe/Surtsey3/'
directorydatap = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directoryfigure = '/home/zlabe/Desktop/RegionalMask/'
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
#sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
#sitp = readPIOMAS(directorydatap,0.15)
#lons2,lats2 = np.meshgrid(lons,lats)
#
#sitall = np.append(sith,sitf,axis=1)

### Slice regions
region = 'Greenland'

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

## See the region for trend calculation
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = fig.add_subplot(111)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,resolution='l',round =True)
m.drawcoastlines(color = 'k',linewidth=0.2)
m.drawcountries(color='k',linewidth=0.2)
m.drawlsmask(land_color='darkgrey',ocean_color='azure')
m.drawmapboundary(color='white')
cs = m.contourf(lons2n,lats2n,lats2n,latlon=True,colors='r')
m.fillcontinents(color='darkgrey')
fig.suptitle(r'\textbf{%s Mask}' % region)
plt.savefig(directoryfigure + 'mask_%s.png' % region,dpi=300)

### Calculate decadal trends
def monRegress(sitq,months,ensemble):
    slopesit = np.zeros((sitq.shape[1],sitq.shape[2],sitq.shape[3]))
    for mo in xrange(sitq.shape[1]):
        sit = sitq[:,mo,:,:]
        for i in xrange(0,sit.shape[1]):
            for j in xrange(0,sit.shape[2]):
                varyy = np.ravel(sit[:,i,j])
                varxx = np.arange(varyy.shape[0])
                mask = np.isfinite(varxx) & np.isfinite(varyy)
                
                varyymean = np.nanmean(varyy)
                if np.isfinite(varyymean):
                    slopesit[mo,i,j],intercept,r,p_value,std_err = sts.stats.linregress(varxx[mask],
                                                                      varyy[mask])
                else:
                    slopesit[mo,i,j] = np.nan  
        print 'Completed: Month %s done!' % (months[mo])
    print 'Completed: Calculated regression!'
    
    slopesit = slopesit*10. # decadal trend
    return slopesit

### Calculate gridded decadal trends
yearq = np.where((years >= 1979) & (years <= 2015))[0]

sittrendhq = np.empty((sitmmh.shape[0],sitmmh.shape[2],sitmmh.shape[3],sitmmh.shape[4]))
sittrendfq = np.empty((sitmmf.shape[0],sitmmf.shape[2],sitmmf.shape[3],sitmmf.shape[4]))
sittrendpq = np.empty((sitmmf.shape[0],sitmmf.shape[2],sitmmf.shape[3],sitmmf.shape[4]))
for i in xrange(sitmmh.shape[0]):
#    sittrendhq[i] = monRegress(sitmmh[i,:,:,:,:],months,ensemble)
    sittrendfq[i] = monRegress(sitmmf[i,:,:,:,:],months,ensemble)
#    sittrendpq[i] = monRegress(sitall[i,yearq,:,:,:],months,ensemble)
    
sittrendPio = monRegress(sitmmp,months,ensemble)

### Select trends
#trendh = sittrendhq
trendf = sittrendfq
#trendp = sittrendpq
trendpio = sittrendPio

### Slice seasons
#trendh_w = np.nanmean(trendh[:,0:3,:,:],axis=1)
#trendh_sp = np.nanmean(trendh[:,3:6,:,:],axis=1)
#trendh_su = np.nanmean(trendh[:,6:9,:,:],axis=1)
#trendh_f = np.nanmean(trendh[:,9:12,:,:],axis=1)

trendf_w = np.nanmean(trendf[:,0:3,:,:],axis=1)
trendf_sp = np.nanmean(trendf[:,3:6,:,:],axis=1)
trendf_su = np.nanmean(trendf[:,6:9,:,:],axis=1)
trendf_f = np.nanmean(trendf[:,9:12,:,:],axis=1)

#trendp_w = np.nanmean(trendp[:,0:3,:,:],axis=1)
#trendp_sp = np.nanmean(trendp[:,3:6,:,:],axis=1)
#trendp_su = np.nanmean(trendp[:,6:9,:,:],axis=1)
#trendp_f = np.nanmean(trendp[:,9:12,:,:],axis=1)

trendpio_w = np.nanmean(trendpio[0:3,:,:],axis=0)
trendpio_sp = np.nanmean(trendpio[3:6,:,:],axis=0)
trendpio_su = np.nanmean(trendpio[6:9,:,:],axis=0)
trendpio_f = np.nanmean(trendpio[9:12,:,:],axis=0)

def weightThick(var,lats,types):
    """
    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
    """
    
    if types == 'lens':
        sityr = np.empty((var.shape[0]))
        for ens in xrange(var.shape[0]):
            varq = var[ens,:,:]
            mask = np.isfinite(varq) & np.isfinite(lats)
            varmask = varq[mask]
            areamask = np.cos(np.deg2rad(lats[mask]))
            sityr[ens] = np.nansum(varmask*areamask)/np.sum(areamask)
            
            print 'Completed: Weighting per ensemble #%s!' % ensemble[ens]
    
    elif types == 'piomas':
        varq = var[:,:]
        mask = np.isfinite(varq) & np.isfinite(lats)
        varmask = varq[mask]
        areamask = np.cos(np.deg2rad(lats[mask]))
        sityr = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr
    
trendmeanf_w = weightThick(trendf_w,lats2n,'lens')
trendmeanf_sp = weightThick(trendf_sp,lats2n,'lens')
trendmeanf_su = weightThick(trendf_su,lats2n,'lens')
trendmeanf_f = weightThick(trendf_f,lats2n,'lens')

trendmeanpio_w = weightThick(trendpio_w,lats2n,'piomas')
trendmeanpio_sp = weightThick(trendpio_sp,lats2n,'piomas')
trendmeanpio_su = weightThick(trendpio_su,lats2n,'piomas')
trendmeanpio_f = weightThick(trendpio_f,lats2n,'piomas')

ense = np.arange(len(ensemble))

### Trends Figure
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
ax = plt.subplot(141)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.axvline(np.nanmean(trendmeanf_w),color='k',linewidth=2,alpha=0.65)
plt.scatter(trendmeanf_w,ense,s=15,color='teal')
plt.axvline(trendmeanpio_w,color='m',linewidth=1.5)

plt.xticks(np.arange(-1,0.1,0.5),
           map(str,np.arange(-1,0.1,0.5)),fontsize=8)
plt.xlim([-1,0])
plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
plt.ylim([0,40])

plt.text(-0.9,40,r'\textbf{JFM}',fontsize=20,color='darkgrey')
plt.ylabel(r'\textbf{Ensemble Number}')

ax = plt.subplot(142)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.axvline(np.nanmean(trendmeanf_sp),color='k',linewidth=2,alpha=0.65)
plt.scatter(trendmeanf_sp,ense,s=15,color='teal')
plt.axvline(trendmeanpio_sp,color='m',linewidth=1.5)

plt.xticks(np.arange(-1,0.1,0.5),
           map(str,np.arange(-1,0.1,0.5)),fontsize=8)
plt.xlim([-1,0])
plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
plt.ylim([0,40])

plt.text(-0.9,40,r'\textbf{AMJ}',fontsize=20,color='darkgrey')

ax = plt.subplot(143)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.axvline(np.nanmean(trendmeanf_su),color='k',linewidth=2,alpha=0.65)
plt.scatter(trendmeanf_su,ense,s=15,color='teal')
plt.axvline(trendmeanpio_su,color='m',linewidth=1.5)

plt.xticks(np.arange(-1,0.1,0.5),
           map(str,np.arange(-1,0.1,0.5)),fontsize=8)
plt.xlim([-1,0])
plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
plt.ylim([0,40])

plt.text(-0.9,40,r'\textbf{JAS}',fontsize=20,color='darkgrey')

ax = plt.subplot(144)
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

plt.axvline(np.nanmean(trendmeanf_f),color='k',linewidth=2,alpha=0.65)
plt.scatter(trendmeanf_f,ense,s=15,color='teal')
plt.axvline(trendmeanpio_f,color='m',linewidth=1.5)

plt.xticks(np.arange(-1,0.1,0.5),
           map(str,np.arange(-1,0.1,0.5)),fontsize=8)
plt.xlim([-1,0])
plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
plt.ylim([0,40])

plt.text(-0.9,40,r'\textbf{OND}',fontsize=20,color='darkgrey')

ax.text(-3.3,-6,r'\textbf{LENS SIT( m decade$^{-1}$ )}')

fig.subplots_adjust(wspace=0.3)
plt.savefig(directoryfigure+'future_%s_lens_sittrends.png' % region,dpi=300)