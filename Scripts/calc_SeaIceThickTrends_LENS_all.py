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
        
#def readPIOMAS(directorydata,threshold):
#    files = 'piomas_regrid_sit_LENS_19792015.nc'
#    filename = directorydata + files
#    
#    data = Dataset(filename)
#    sitp = data.variables['sit'][:,:,156:180,:] # lats > 65
#    data.close()
#    
#    ### Mask out threshold values
#    if threshold == 'None':
#        sitp[np.where(sitp < 0)] = np.nan
#        sitp[np.where(sitp > 12)] = np.nan
#    else:
#        sitp[np.where(sitp < threshold)] = np.nan
#        sitp[np.where(sitp < 0)] = np.nan
#        sitp[np.where(sitp > 12)] = np.nan
#    
#    print 'Completed: Read PIOMAS SIT!'
#    return sitp
#    
#### Call functions   
#sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
#sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
#sitp = readPIOMAS(directorydatap,0.15)
#lons,lats = np.meshgrid(lons,lats)
#
#sitall = np.append(sith,sitf,axis=1)
#
#### Calculate decadal trends
#def monRegress(sitq,months):
#    slopesit = np.zeros((sitq.shape[1],sitq.shape[2],sitq.shape[3]))
#    for mo in xrange(sitq.shape[1]):
#        sit = sitq[:,mo,:,:]
#        for i in xrange(0,sit.shape[1]):
#            for j in xrange(0,sit.shape[2]):
#                varyy = np.ravel(sit[:,i,j])
#                varxx = np.arange(varyy.shape[0])
#                mask = np.isfinite(varxx) & np.isfinite(varyy)
#                
#                varyymean = np.nanmean(varyy)
#                if np.isfinite(varyymean):
#                    slopesit[mo,i,j],intercept,r,p_value,std_err = sts.stats.linregress(varxx[mask],
#                                                                      varyy[mask])
#                else:
#                    slopesit[mo,i,j] = np.nan  
#        print 'Completed: Month %s done!' % (months[mo])
#    print 'Completed: Calculated regression!'
#    
#    slopesit = slopesit*10. # decadal trend
#    return slopesit
#
#### Calculate gridded decadal trends
#yearq = np.where((years >= 1979) & (years <= 2015))[0]
#
#sittrendhq = np.empty((sith.shape[0],sith.shape[2],sith.shape[3],sith.shape[4]))
#sittrendfq = np.empty((sitf.shape[0],sitf.shape[2],sitf.shape[3],sitf.shape[4]))
#sittrendpq = np.empty((sitf.shape[0],sitf.shape[2],sitf.shape[3],sitf.shape[4]))
#for i in xrange(sitall.shape[0]):
#    sittrendhq[i] = monRegress(sith[i,:,:,:,:],months)
#    sittrendfq[i] = monRegress(sitf[i,:,:,:,:],months)
#    sittrendpq[i] = monRegress(sitall[i,yearq,:,:,:],months)
#    
#sittrendPio = monRegress(sitp,months)

#def netcdfPiomas(lats,lons,var,directory):
#    print '\n>>> Using netcdf4LENS function!'
#    
#    name = 'piomas_sittrend_19792015.nc'
#    filename = directory + name
#    ncfile = Dataset(filename,'w',format='NETCDF4')
#    ncfile.description = 'piomas decadal trend sit interpolated on 1x1 grid' 
#    
#    ### Dimensions
##    ncfile.createDimension('ensemble',var.shape[0])
#    ncfile.createDimension('months',var.shape[0])
#    ncfile.createDimension('lat',var.shape[1])
#    ncfile.createDimension('lon',var.shape[2])
#    
#    ### Variables
##    ensemble = ncfile.createVariable('ensemble','f4',('ensemble'))
#    months = ncfile.createVariable('months','f4',('months'))
#    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
#    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
#    varns = ncfile.createVariable('trend','f4',('months','lat','lon'))
#    
#    ### Units
#    varns.units = 'meters'
#    ncfile.title = 'LENS sit decadal trend'
#    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
#    ncfile.source = 'University of Washington, Polar Science Center'
#    ncfile.references = 'Zhang and Rothrock [2003]'
#    
#    ### Data
##    ensemble[:] = list(xrange(var.shape[0]))
#    months[:] = list(xrange(var.shape[0]))
#    latitude[:] = lats
#    longitude[:] = lons
#    varns[:] = var
#    
#    ncfile.close()
#    print '*Completed: Created netCDF4 File!'

#netcdfPiomas(lats,lons,sittrendhq,directorydata2)
#netcdfPiomas(lats,lons,sittrendfq,directorydata2)
#netcdfPiomas(lats,lons,sittrendpq,directorydata2)
#netcdfPiomas(lats,lons,sittrendPio,directorydata2)

###########################################################################
###########################################################################
###########################################################################
### Read in trends
data = Dataset(directorydata2 + 'lens_sittrend_19202005.nc')
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
trendh = data.variables['trend'][:]
data.close()

data = Dataset(directorydata2 + 'lens_sittrend_20062080.nc')
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
trendf = data.variables['trend'][:]
data.close()

data = Dataset(directorydata2 + 'lens_sittrend_19792015.nc')
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
trendp = data.variables['trend'][:]
data.close()

data = Dataset(directorydata2 + 'piomas_sittrend_19792015.nc')
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
trendpio = data.variables['trend'][:]
data.close()

### Slice seasons
trendh_w = np.nanmean(trendh[:,0:3,:,:],axis=1)
trendh_sp = np.nanmean(trendh[:,3:6,:,:],axis=1)
trendh_su = np.nanmean(trendh[:,6:9,:,:],axis=1)
trendh_f = np.nanmean(trendh[:,9:12,:,:],axis=1)

trendf_w = np.nanmean(trendf[:,0:3,:,:],axis=1)
trendf_sp = np.nanmean(trendf[:,3:6,:,:],axis=1)
trendf_su = np.nanmean(trendf[:,6:9,:,:],axis=1)
trendf_f = np.nanmean(trendf[:,9:12,:,:],axis=1)

trendp_w = np.nanmean(trendp[:,0:3,:,:],axis=1)
trendp_sp = np.nanmean(trendp[:,3:6,:,:],axis=1)
trendp_su = np.nanmean(trendp[:,6:9,:,:],axis=1)
trendp_f = np.nanmean(trendp[:,9:12,:,:],axis=1)

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
    
trendmeanf_w = weightThick(trendf_w,lats,'lens')
trendmeanf_sp = weightThick(trendf_sp,lats,'lens')
trendmeanf_su = weightThick(trendf_su,lats,'lens')
trendmeanf_f = weightThick(trendf_f,lats,'lens')

trendmeanpio_w = weightThick(trendpio_w,lats,'piomas')
trendmeanpio_sp = weightThick(trendpio_sp,lats,'piomas')
trendmeanpio_su = weightThick(trendpio_su,lats,'piomas')
trendmeanpio_f = weightThick(trendpio_f,lats,'piomas')

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
        
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

#fig = plt.figure()
#ax = plt.subplot(141)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.axvline(np.nanmean(trendmeanf_w),color='k',linewidth=2,alpha=0.65)
#plt.scatter(trendmeanf_w,ense,s=15,color='teal')
#plt.axvline(trendmeanpio_w,color='m',linewidth=1.5)
#
#plt.xticks(np.arange(-0.5,0.1,0.25),
#           map(str,np.arange(-0.5,0.1,0.25)),fontsize=8)
#plt.xlim([-0.5,0])
#plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
#plt.ylim([0,40])
#
#plt.text(-0.4,40,r'\textbf{JFM}',fontsize=20,color='darkgrey')
#plt.ylabel(r'\textbf{Ensemble Number}')
#
#ax = plt.subplot(142)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.axvline(np.nanmean(trendmeanf_sp),color='k',linewidth=2,alpha=0.65)
#plt.scatter(trendmeanf_sp,ense,s=15,color='teal')
#plt.axvline(trendmeanpio_sp,color='m',linewidth=1.5)
#
#plt.xticks(np.arange(-0.5,0.1,0.25),
#           map(str,np.arange(-0.5,0.1,0.25)),fontsize=8)
#plt.xlim([-0.5,0])
#plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
#plt.ylim([0,40])
#
#plt.text(-0.4,40,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
#
#ax = plt.subplot(143)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.axvline(np.nanmean(trendmeanf_su),color='k',linewidth=2,alpha=0.65)
#plt.scatter(trendmeanf_su,ense,s=15,color='teal')
#plt.axvline(trendmeanpio_su,color='m',linewidth=1.5)
#
#plt.xticks(np.arange(-0.5,0.1,0.25),
#           map(str,np.arange(-0.5,0.1,0.25)),fontsize=8)
#plt.xlim([-0.5,0])
#plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
#plt.ylim([0,40])
#
#plt.text(-0.4,40,r'\textbf{JAS}',fontsize=20,color='darkgrey')
#
#ax = plt.subplot(144)
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.axvline(np.nanmean(trendmeanf_f),color='k',linewidth=2,alpha=0.65)
#plt.scatter(trendmeanf_f,ense,s=15,color='teal')
#plt.axvline(trendmeanpio_f,color='m',linewidth=1.5)
#
#plt.xticks(np.arange(-0.5,0.1,0.25),
#           map(str,np.arange(-0.5,0.1,0.25)),fontsize=8)
#plt.xlim([-0.5,0])
#plt.yticks(np.arange(0,45,5),map(str,np.arange(0,45,5)),fontsize=8)
#plt.ylim([0,40])
#
#plt.text(-0.4,40,r'\textbf{OND}',fontsize=20,color='darkgrey')
#
#ax.text(-1.6,-6,r'\textbf{LENS SIT( m decade$^{-1}$ )}')
#
#fig.subplots_adjust(wspace=0.3)
#plt.savefig(directoryfigure+'future_lens_sittrends.png',dpi=300)

### Figures
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

timex = np.arange(-1,1.5,0.5)
timey = np.arange(-1,1.5,0.5)

varx_w = trendpio_w
vary_w = np.nanmean(trendp_w,axis=0)      
mask = np.isfinite(varx_w) & np.isfinite(vary_w)
slope_w, intercept_w, r_value_w, p_value, std_err = sts.linregress(varx_w[mask], vary_w[mask])
line_w = slope_w*timex + intercept_w

varx_sp = trendpio_sp
vary_sp = np.nanmean(trendp_sp,axis=0)      
mask = np.isfinite(varx_sp) & np.isfinite(vary_sp)
slope_sp, intercept_sp, r_value_sp, p_value, std_err = sts.linregress(varx_sp[mask], vary_sp[mask])
line_sp = slope_sp*timex + intercept_sp

varx_su = trendpio_su
vary_su = np.nanmean(trendp_su,axis=0)      
mask = np.isfinite(varx_su) & np.isfinite(vary_su)
slope_su, intercept_su, r_value_su, p_value, std_err = sts.linregress(varx_su[mask], vary_su[mask])
line_su = slope_su*timex + intercept_su

varx_f = trendpio_w
vary_f = np.nanmean(trendp_f,axis=0)      
mask = np.isfinite(varx_f) & np.isfinite(vary_f)
slope_f, intercept_f, r_value_f, p_value, std_err = sts.linregress(varx_f[mask], vary_f[mask])
line_f = slope_f*timex + intercept_f
        
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(221)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

ax.plot(timex,timey,color='k',linewidth=3,zorder=1)
for i in xrange(trendp_w.shape[0]):
    plt.scatter(trendpio_w,trendp_w[i],color='darkgrey',s=0.1,alpha=0.5,
                zorder=2)
plt.scatter(trendpio_w,np.nanmean(trendp_w,axis=0),color='teal',s=1,
            zorder=3,alpha=0.5)
plt.plot(timex,line_w,color='r',zorder=4)

plt.xticks(np.arange(-1,1.5,0.5),
           map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.xlim([-1,1])
plt.yticks(np.arange(-1,1.5,0.5),map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.ylim([-1,1])

plt.ylabel(r'\textbf{LENS SIT( m decade$^{-1}$ )}')
plt.text(-1.05,0.9,r'\textbf{JFM}',fontsize=20,color='darkgrey')
plt.text(0.5,-1,r'LENS Mean',fontsize=7,color='teal')
plt.text(-1.05,0.7,r'R$^{2}$=%s' % round(r_value_w**2,2),fontsize=8,color='r')

ax = plt.subplot(222)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

ax.plot(timex,timey,color='k',linewidth=3,zorder=1)
for i in xrange(trendp_w.shape[0]):
    plt.scatter(trendpio_sp,trendp_sp[i],color='darkgrey',s=0.1,alpha=0.5,
                zorder=2)
plt.scatter(trendpio_sp,np.nanmean(trendp_sp,axis=0),color='teal',s=1,
            zorder=3,alpha=0.5)
plt.plot(timex,line_sp,color='r',zorder=4)

plt.xticks(np.arange(-1,1.5,0.5),
           map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.xlim([-1,1])
plt.yticks(np.arange(-1,1.5,0.5),map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.ylim([-1,1])

plt.text(-1.05,0.9,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
plt.text(0.5,-1,r'LENS Mean',fontsize=7,color='teal')
plt.text(-1.05,0.7,r'R$^{2}$=%s' % round(r_value_sp**2,2),fontsize=8,color='r')

ax = plt.subplot(223)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

ax.plot(timex,timey,color='k',linewidth=3,zorder=1)
for i in xrange(trendp_w.shape[0]):
    plt.scatter(trendpio_su,trendp_su[i],color='darkgrey',s=0.1,alpha=0.5,
                zorder=2)
plt.scatter(trendpio_su,np.nanmean(trendp_su,axis=0),color='teal',s=1,
            zorder=3,alpha=0.5)
plt.plot(timex,line_su,color='r',zorder=4)

plt.xticks(np.arange(-1,1.5,0.5),
           map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.xlim([-1,1])
plt.yticks(np.arange(-1,1.5,0.5),map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.ylim([-1,1])

plt.text(-1.05,0.9,r'\textbf{JAS}',fontsize=20,color='darkgrey')
plt.text(0.5,-1,r'LENS Mean',fontsize=7,color='teal')
plt.text(-1.05,0.7,r'R$^{2}$=%s' % round(r_value_su**2,2),fontsize=8,color='r')

plt.xlabel(r'\textbf{PIOMAS SIT( m decade$^{-1}$ )}')
plt.ylabel(r'\textbf{LENS SIT( m decade$^{-1}$ )}')

ax = plt.subplot(224)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

ax.plot(timex,timey,color='k',linewidth=3,zorder=1)
for i in xrange(trendp_f.shape[0]):
    plt.scatter(trendpio_f,trendp_f[i],color='darkgrey',s=0.1,alpha=0.5,
                zorder=2)
plt.scatter(trendpio_f,np.nanmean(trendp_f,axis=0),color='teal',s=1,
            zorder=3,alpha=0.5)
plt.plot(timex,line_f,color='r',zorder=4)

plt.xticks(np.arange(-1,1.5,0.5),
           map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.xlim([-1,1])
plt.yticks(np.arange(-1,1.5,0.5),map(str,np.arange(-1,1.5,0.5)),fontsize=8)
plt.ylim([-1,1])

plt.text(-1.05,0.9,r'\textbf{OND}',fontsize=20,color='darkgrey')
plt.text(0.5,-1,r'LENS Mean',fontsize=7,color='teal')
plt.text(-1.05,0.7,r'R$^{2}$=%s' % round(r_value_f**2,2),fontsize=8,color='r')

plt.xlabel(r'\textbf{PIOMAS SIT( m decade$^{-1}$ )}')

fig.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'LENSPIOMAS_trends_scatter.png',dpi=300)