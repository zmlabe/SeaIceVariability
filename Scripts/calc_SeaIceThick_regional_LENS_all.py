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
#    sitmp = sitp[:,:,latq,:]
    
    # Slice lons
    sitmmh = sitmh[:,:,:,:,lonq]
    sitmmf = sitmf[:,:,:,:,lonq]
#    sitmmp = sitmp[:,:,:,lonq]
    
    sitall = np.append(sitmmh,sitmmf,axis=1)
    
    return sitall,lats2n,lons2n
    
#### Regions sliced   
#sitallg,lats2ng,lons2ng = RegionalSlice(sith,sitf,'Greenland',lats,lons)
#sitallb,lats2nb,lons2nb = RegionalSlice(sith,sitf,'BeaufortSea',lats,lons)
#sitalle,lats2ne,lons2ne = RegionalSlice(sith,sitf,'EastSiberianSea',lats,lons)
#sitalll,lats2nl,lons2nl = RegionalSlice(sith,sitf,'LaptevSea',lats,lons)
#sitallkb,lats2nkb,lons2nkb = RegionalSlice(sith,sitf,'KB',lats,lons)
#sitallcab,lats2ncab,lons2ncab = RegionalSlice(sith,sitf,'CAB',lats,lons)

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
#sityrg = weightThick(sitallg,lats2ng,'lens')    
#sityrb = weightThick(sitallb,lats2nb,'lens') 
#sityre = weightThick(sitalle,lats2ne,'lens') 
#sityrl = weightThick(sitalll,lats2nl,'lens') 
#sityrkb = weightThick(sitallkb,lats2nkb,'lens') 
#sityrcab = weightThick(sitallcab,lats2ncab,'lens')     
#
#### Select month (mar/sep)
#monthq = 8
#if monthq < 12:
#    sityrg = sityrg[:,:,monthq]
#    sityrb = sityrb[:,:,monthq]
#    sityre = sityre[:,:,monthq]
#    sityrl = sityrl[:,:,monthq]
#    sityrkb = sityrkb[:,:,monthq]
#    sityrcab = sityrcab[:,:,monthq]
#
#### Calculate mean per region
#sitmeang = np.nanmean(sityrg,axis=0)  
#sitmeanb = np.nanmean(sityrb,axis=0)  
#sitmeane = np.nanmean(sityre,axis=0)  
#sitmeanl = np.nanmean(sityrl,axis=0)   
#sitmeankb = np.nanmean(sityrkb,axis=0)  
#sitmeancab = np.nanmean(sityrcab,axis=0)   
#        
##### Plot Figure
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#### Adjust axes in time series plots 
#def adjust_spines(ax, spines):
#    for loc, spine in ax.spines.items():
#        if loc in spines:
#            spine.set_position(('outward', 10))
#        else:
#            spine.set_color('none')  
#    if 'left' in spines:
#        ax.yaxis.set_ticks_position('left')
#    else:
#        ax.yaxis.set_ticks([])
#
#    if 'bottom' in spines:
#        ax.xaxis.set_ticks_position('bottom')
#    else:
#        ax.xaxis.set_ticks([]) 
#        
#fig = plt.figure()
#ax = plt.subplot(111)
#
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#         
#plt.plot(sitmeang,color='steelblue',linewidth=2,linestyle='-',
#         zorder=3,label=r'Greenland Sea')    
#plt.plot(sitmeanb,color='darkgreen',linewidth=2,linestyle='-',
#         zorder=4,label=r'Beaufort Sea')   
#plt.plot(sitmeane,color='darkorange',linewidth=2,linestyle='-',
#         zorder=5,label=r'East Siberian Sea')
#plt.plot(sitmeanl,color='darkblue',linewidth=2,linestyle='-',
#         zorder=6,label=r'Laptev Sea')    
#plt.plot(sitmeankb,color='m',linewidth=2,linestyle='-',
#         zorder=7,label=r'Kara-Barents Seas')   
#plt.plot(sitmeancab,color='yellowgreen',linewidth=2,linestyle='-',
#         zorder=8,label=r'Central Arctic Basin')
#         
#plt.axvline(85,linestyle='--',linewidth=2,color='k')
#
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.yticks(np.arange(0,6,1),map(str,np.arange(0,6,1))) 
#plt.xlim([0,160])
#plt.ylim([0,5])
#
#plt.ylabel(r'\textbf{Sea Ice Thickness (m)}')
#
#plt.legend(shadow=False,fontsize=7,loc='upper right',
#           fancybox=True,frameon=False,ncol=1)
#           
#plt.savefig(directoryfigure+'lens_yr_regional_sep_all.png',dpi=300)
#
############################################################################
############################################################################
############################################################################
#
#fig = plt.figure(figsize=(5,9))
#
#ax = plt.subplot(612)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#for i in xrange(sityrg.shape[0]):
#    plt.plot(sityrg[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeang,color='steelblue',linewidth=2,linestyle='-',
#         zorder=3,label=r'Greenland Sea')  
#         
#### Below 0.5 m
#thresg = np.where(sitmeang <= 0.5)[0]
#if thresg.shape[0] > 0:
#    plt.axvline(thresg[0],linewidth=2,color='k')
#         
##xlabels = map(str,np.arange(1979,2017,5))
##plt.xticks(np.arange(0,39,5),xlabels)
#plt.xlim([0,160])
#
#plt.text(0.5,0,r'\textbf{GREENLAND}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#         
############################################################################         
#         
#ax = plt.subplot(613)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#for i in xrange(sityrb.shape[0]):
#    plt.plot(sityrb[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeanb,color='darkgreen',linewidth=2,linestyle='-',
#         zorder=4,label=r'Beaufort Sea')
#         
#### Below 0.5 m
#thresb = np.where(sitmeanb <= 0.5)[0]
#if thresb.shape[0] > 0:
#    plt.axvline(thresb[0],linewidth=2,color='k')
#
##xlabels = map(str,np.arange(1979,2017,5))
##plt.xticks(np.arange(0,39,5),xlabels)
#plt.xlim([0.5,160])
#
#plt.text(0,0,r'\textbf{BEAUFORT-CHUKCHI SEAS}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#         
############################################################################          
#         
#ax = plt.subplot(614)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#for i in xrange(sityre.shape[0]):
#    plt.plot(sityre[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeane,color='darkorange',linewidth=2,linestyle='-',
#         zorder=5,label=r'East Siberian Sea')
#         
#### Below 0.5 m
#threse = np.where(sitmeane <= 0.5)[0]
#if threse.shape[0] > 0:
#    plt.axvline(threse[0],linewidth=2,color='k')
#
##xlabels = map(str,np.arange(1979,2017,5))
##plt.xticks(np.arange(0,39,5),xlabels)
#plt.xlim([0,160])
#
#plt.text(0.5,0,r'\textbf{EAST SIBERIAN SEA}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#
############################################################################ 
#
#ax = plt.subplot(615)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#for i in xrange(sityrl.shape[0]):
#    plt.plot(sityrl[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeanl,color='darkblue',linewidth=2,linestyle='-',
#         zorder=6,label=r'Laptev Sea')   
#         
#### Below 0.5 m
#thresl = np.where(sitmeanl <= 0.5)[0]
#if thresl.shape[0] > 0:
#    plt.axvline(thresl[0],linewidth=2,color='k')
#
##xlabels = map(str,np.arange(1979,2017,5))
##plt.xticks(np.arange(0,39,5),xlabels)
#plt.xlim([0.5,160])
#
#plt.text(0,0,r'\textbf{LAPTEV SEA}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#
############################################################################ 
#
#ax = plt.subplot(616)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=4,width=1.5,which='major',color='darkgrey')
#
#for i in xrange(sityrkb.shape[0]):
#    plt.plot(sityrkb[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeankb,color='m',linewidth=2,linestyle='-',
#         zorder=7,label=r'Kara-Barents Seas')   
#
#### Below 0.5 m
#threskb = np.where(sitmeankb <= 0.5)[0]
#if threskb.shape[0] > 0:
#    plt.axvline(threskb[0],linewidth=2,color='k')
#
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.xlim([0,160])
#
#plt.text(0.5,0,r'\textbf{BARENTS-KARA SEAS}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#
############################################################################ 
#
#ax = plt.subplot(611)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#for i in xrange(sityrcab.shape[0]):
#    plt.plot(sityrcab[i,:],color='dimgrey',alpha=0.55,linewidth=0.2,
#         zorder=1)
#plt.plot(sitmeancab,color='r',linewidth=2,linestyle='-',
#         zorder=8,label=r'Central Arctic Basin')
#
#### Below 0.5 m
#threscab = np.where(sitmeancab <= 0.5)[0]
#if threscab.shape[0] > 0:
#    plt.axvline(threscab[0],linewidth=2,color='k')
#
#plt.text(0.5,0,r'\textbf{CENTRAL ARCTIC BASIN}',fontsize=11,color='darkgrey')
#plt.text(-22,-9.4,r'\textbf{Sea Ice Thickness (m)}',fontsize=15,
#                          color='k',rotation=90)
#
#plt.yticks(np.arange(0,7,2),map(str,np.arange(0,7,2))) 
#plt.ylim([0,6])
#
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.savefig(directoryfigure + 'lens_yr_regional_sep_all_2.png',dpi=300)
#
############################################################################
############################################################################
############################################################################

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

latmin = 76
latmax = 85
lonmin = 225
lonmax = 370
lonmin2 = 0
lonmax2 = 20
latq = np.where((lats >= latmin) & (lats <= latmax))[0]
lats1ng = lats[latq]
lonq1 = np.where((lons >= lonmin) & (lons <= lonmax))[0]
lonq2 = np.where((lons >= lonmin2) & (lons <= lonmax2))[0]
lonq = np.append(lonq2,lonq1,axis=0)
lons1ng = lons[lonq]
lons2ng,lats2ng = np.meshgrid(lons1ng,lats1ng)

var, lons_cyclic = addcyclic(lats2ng, lons1ng)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats1ng)
x, y = m(lon2d, lat2d)    

parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
m.drawparallels(parallels,labels=[False,False,False,False],linewidth=0,color='k')
m.drawmeridians(meridians,labels=[True,True,False,False],linewidth=0,fontsize=4,color='k')

cs = m.contourf(x,y,var,colors='steelblue')
cs = m.contourf(lons2nb,lats2nb,lats2nb,latlon=True,colors='darkgreen')
cs = m.contourf(lons2ne,lats2ne,lats2ne,latlon=True,colors='darkorange')
cs = m.contourf(lons2nl,lats2nl,lats2nl,latlon=True,colors='darkblue')
cs = m.contourf(lons2nkb,lats2nkb,lats2nkb,latlon=True,colors='m')
cs = m.contourf(lons2ncab,lats2ncab,lats2ncab,latlon=True,colors='salmon')
m.fillcontinents(color='darkgrey')

fig.text(0.58,0.63,r'\textbf{B-K}',color='k',fontsize=20)
fig.text(0.462,0.48,r'\textbf{CAB}',color='k',fontsize=20)
fig.text(0.35,0.53,r'\textbf{ESS}',color='k',fontsize=20)
fig.text(0.35,0.4,r'\textbf{B-C}',color='k',fontsize=20)
fig.text(0.425,0.63,r'\textbf{LV}',color='k',fontsize=20)
fig.text(0.515,0.365,r'\textbf{GD}',color='k',fontsize=20)

plt.savefig(directoryfigure + 'masking_regions.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
#
#fig = plt.figure(figsize=(5,9))
#
#ax = plt.subplot(611)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#varig = np.nanmax(sityrg,axis=0) - np.nanmin(sityrg,axis=0)
#
#plt.plot(varig,color='steelblue',linewidth=2,linestyle='-',
#         zorder=3,label=r'Greenland Sea')   
#         
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{GREENLAND}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#         
############################################################################         
#         
#ax = plt.subplot(612)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#varib = np.nanmax(sityrb,axis=0) - np.nanmin(sityrb,axis=0)
#
#plt.plot(varib,color='darkgreen',linewidth=2,linestyle='-',
#         zorder=4,label=r'Beaufort Sea')
#
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{BEAUFORT-CHUKCHI SEAS}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#         
############################################################################          
#         
#ax = plt.subplot(613)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#varie = np.nanmax(sityre,axis=0) - np.nanmin(sityre,axis=0)
#         
#plt.plot(varie,color='darkorange',linewidth=2,linestyle='-',
#         zorder=5,label=r'East Siberian Sea')
#
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{EAST SIBERIAN SEA}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#
############################################################################ 
#
#ax = plt.subplot(614)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#varil = np.nanmax(sityrl,axis=0) - np.nanmin(sityrl,axis=0)
#         
#plt.plot(varil,color='darkblue',linewidth=2,linestyle='-',
#         zorder=6,label=r'Laptev Sea')   
#
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{LAPTEV SEA}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#
############################################################################ 
#
#ax = plt.subplot(615)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('none')
#ax.tick_params('y',length=4,width=1.5,which='major',color='darkgrey')
#ax.tick_params('x',length=0,width=0,which='major',color='none')
#ax.xaxis.set_visible(False)
#
#varikb = np.nanmax(sityrkb,axis=0) - np.nanmin(sityrkb,axis=0)
#         
#plt.plot(varikb,color='m',linewidth=2,linestyle='-',
#         zorder=7,label=r'Kara-Barents Seas')   
#
##plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
##plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{BARENTS-KARA SEAS}',fontsize=11,color='darkgrey')
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#
############################################################################ 
#
#ax = plt.subplot(616)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#varicab = np.nanmax(sityrcab,axis=0) - np.nanmin(sityrcab,axis=0)
#
#plt.plot(varicab,color='r',linewidth=2,linestyle='-',
#         zorder=8,label=r'Central Arctic Basin')
#         
#plt.xticks(np.arange(0,181,20),np.arange(1920,2101,20))
#plt.xlim([0,160])
#
#plt.text(0,3.5,r'\textbf{CENTRAL ARCTIC BASIN}',fontsize=11,color='darkgrey')
#plt.text(-22,16.3,r'\textbf{Difference (m)}',fontsize=15,
#                          color='k',rotation=90)
#
#plt.yticks(np.arange(0,5,2),map(str,np.arange(0,5,2))) 
#plt.ylim([0,4])
#
#plt.savefig(directoryfigure + 'SIT_regions_diff_maxmin_sep.png',dpi=300)
#
#############################################################################
#############################################################################
#############################################################################
