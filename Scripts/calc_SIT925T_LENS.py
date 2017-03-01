"""
Script calculates comparisons of T2M with SIT from LENS. Note the first 
part regrids on 1x1 degree grid for comparison with LENS SIT
 
Notes
-----
    Author : Zachary Labe
    Date   : 22 February 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_var_LENS as LV
import read_SeaIceThick_LENS as lens
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm
from scipy.interpolate import griddata as g
from netCDF4 import Dataset

### Define directories
directorydataT = '/home/zlabe/Surtsey3/CESM_large_ensemble/'
directorydataSIT = '/home/zlabe/Surtsey3/' 
directorydataN = '/home/zlabe/Surtsey/LENS/T2M/'
directoryfigure = '/home/zlabe/Desktop/'
directorydata2 = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot LENS AD - %s----' % titletime 

ense = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 1920
year2 = 2080
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

yearslens = np.arange(1920,2080+1,1)          
yearsclimo = np.arange(1981,2010+1,1)
          
### Read in functions
##tas,lats1,lons1 = LV.readLENSEnsemble(directorydataT,'T2M') 
#sith,lats2,lons2 = lens.readLENSEnsemble(directorydataSIT,0.15,'historical')
#sitf,lats2,lons2 = lens.readLENSEnsemble(directorydataSIT,0.15,'rcp85')
#
#### Combine SIT periods
#sitall = np.append(sith,sitf,axis=1)
#          
##### 2D lat/lon arrays          
##lons2,lats2 = np.meshgrid(lons2,lats2)
##lons1,lats1 = np.meshgrid(lons1,lats1)
#          
############################################################################  
############################################################################
############################################################################
##### Regrid
##def regrid(lat1,lon1,lat2,lon2,var,years):
##    """
##    Interpolated on selected grid. 
##    [year,month,lat,lon]
##    """
##    
##    varn_re = np.reshape(var,(var.shape[0],var.shape[1],(11*144)))   
##    
##    varn = np.empty((var.shape[0],var.shape[1],lat2.shape[0],lon2.shape[1]))
##    
##    print 'Completed: Start regridding process:'
##    
##    for i in xrange(varn.shape[0]):
##        for j in xrange(varn.shape[1]):
##            z = g((np.ravel(lat1),np.ravel(lon1)),varn_re[i,j,:],(lat2,lon2),method='linear')
##            varn[i,j,:,:] = z
##        print 'Completed: Year %s Regridding---' % (years[i])
##    return varn
#    
##tasn = np.empty(sit.shape)
##for i in xrange(len(ense)):
##    tasn[i,:,:,:,:] = regrid(lats1,lons1,lats2,lons2,tas[i],years)
##    print 'Completed: Regridded Ensemble #%s!' % ense[i]
#
##def netcdfPiomas(lats,lons,var,directory):
##    print '\n>>> Using netcdf4LENS function!'
##    
##    name = 'T2M/lens_regrid_T2M_19202080.nc'
##    filename = directory + name
##    ncfile = Dataset(filename,'w',format='NETCDF4')
##    ncfile.description = 'LENS T2M interpolated on 1x1 grid' 
##    
##    ### Dimensions
##    ncfile.createDimension('ensemble',var.shape[0])
##    ncfile.createDimension('years',var.shape[1])
##    ncfile.createDimension('months',var.shape[2])
##    ncfile.createDimension('lat',var.shape[3])
##    ncfile.createDimension('lon',var.shape[4])
##    
##    ### Variables
##    ensemble = ncfile.createVariable('ensemble','f4',('ensemble'))
##    years = ncfile.createVariable('years','f4',('years'))
##    months = ncfile.createVariable('months','f4',('months'))
##    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
##    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
##    varns = ncfile.createVariable('T2M','f4',('ensemble','years','months','lat','lon'))
##    
##    ### Units
##    varns.units = 'Degrees Celsius'
##    ncfile.title = 'LENS T2M'
##    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
##    ncfile.source = 'NCAR LENS'
##    ncfile.references = 'Kay et al. [2013]'
##    
##    ### Data
##    ensemble[:] = list(xrange(var.shape[0]))
##    years[:] = list(xrange(var.shape[1]))
##    months[:] = list(xrange(var.shape[2]))
##    latitude[:] = lats
##    longitude[:] = lons
##    varns[:] = var
##    
##    ncfile.close()
##    print '*Completed: Created netCDF4 File!'
##
##netcdfPiomas(lats2,lons2,tasn,directorydataN)
#    
############################################################################
############################################################################
############################################################################
#### Read T2M
#data = Dataset(directorydataN + 'lens_regrid_T2M_19202080.nc')
#tasall = data.variables['T2M'][:]
#data.close()

##print 'Completed: Read T2M data! \n'
#
#### Calculate ensemble average
##tasmean = np.nanmean(tas,axis=0)
##sitmean = np.nanmean(sit,axis=0)
##
##print 'Completed: Calculated ensemble mean! \n'
##
##def netcdfPiomas(lats,lons,var,directory):
##    print '\n>>> Using netcdf4LENS function!'
##    
##    name = 'lensmean_regrid_sit_19202080.nc'
##    filename = directory + name
##    ncfile = Dataset(filename,'w',format='NETCDF4')
##    ncfile.description = 'LENS mean sit interpolated on 1x1 grid' 
##    
##    ### Dimensions
##    ncfile.createDimension('years',var.shape[0])
##    ncfile.createDimension('months',var.shape[1])
##    ncfile.createDimension('lat',var.shape[2])
##    ncfile.createDimension('lon',var.shape[3])
##    
##    ### Variables
##    years = ncfile.createVariable('years','f4',('years'))
##    months = ncfile.createVariable('months','f4',('months'))
##    latitude = ncfile.createVariable('lat','f4',('lat','lon'))
##    longitude = ncfile.createVariable('lon','f4',('lat','lon'))
##    varns = ncfile.createVariable('sit','f4',('years','months','lat','lon'))
##    
##    ### Units
##    varns.units = 'meters'
##    ncfile.title = 'LENS sit'
##    ncfile.instituion = 'Dept. ESS at University of California, Irvine'
##    ncfile.source = 'University of Washington, Polar Science Center'
##    ncfile.references = 'Zhang and Rothrock [2003]'
##    
##    ### Data
##    years[:] = list(xrange(var.shape[0]))
##    months[:] = list(xrange(var.shape[1]))
##    latitude[:] = lats
##    longitude[:] = lons
##    varns[:] = var
##    
##    ncfile.close()
##    print '*Completed: Created netCDF4 File!'
##
##netcdfPiomas(lats2,lons2,sitmean,directorydata2)
#
############################################################################
############################################################################
############################################################################
#### Read in ensemble means of sit and T2M
#data = Dataset(directorydata2 + 'lensmean_regrid_sit_19202080.nc')
#lats = data.variables['lat'][:]
#lons = data.variables['lon'][:]
#sit = data.variables['sit'][:]
#data.close()
#
#data = Dataset(directorydata2 + 'lensmean_regrid_T2M_19202080.nc')
#lats = data.variables['lat'][:]
#lons = data.variables['lon'][:]
#tas = data.variables['T2M'][:]
#data.close()
#
#print '\nCompleted: Read in LENS means! \n'
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
##yearq = np.where((years >= 1979) & (years <= 2015))[0]
#
##slopesithq = monRegress(sit[yearq,:,:,:],months)  # 86 and 75
##slopesitfq = monRegress(sit[86:,:,:,:],months)
##
##slopetashq = monRegress(tas[yearq,:,:,:],months)
##slopetasfq = monRegress(tas[86:,:,:,:],months)
#
#def weightAve(var,lats):
#    """
#    Area weights sit array 4d [year,month,lat,lon] into [year,month]
#    """
#    sityr = np.empty((var.shape[0],var.shape[1]))
#    for i in xrange(var.shape[0]):
#        for j in xrange(var.shape[1]):
#            varq = var[i,j,:,:]
#            mask = np.isfinite(varq) & np.isfinite(lats)
#            varmask = varq[mask]
#            areamask = np.cos(np.deg2rad(lats[mask]))
#            sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
#     
#    print '\nCompleted: Yearly weighted SIT average!' 
#    return sityr
#
##### Call functions     
#slopesith = weightAve(sit[:86],lats)
#slopesitf = weightAve(sit[86:],lats)
#slopetash = weightAve(tas[:86],lats)
#slopetasf = weightAve(tas[86:],lats)
#
#### Calculate trends for seasons
##sittrendh_w = np.nanmean(slopesith[0:3],axis=0)
##sittrendh_sp = np.nanmean(slopesith[3:6],axis=0)
##sittrendh_su = np.nanmean(slopesith[6:9],axis=0)
##sittrendh_f = np.nanmean(slopesith[9:12],axis=0)
##
##sittrendf_w = np.nanmean(slopesitf[0:3],axis=0)
##sittrendf_sp = np.nanmean(slopesitf[3:6],axis=0)
##sittrendf_su = np.nanmean(slopesitf[6:9],axis=0)
##sittrendf_f = np.nanmean(slopesitf[9:12],axis=0)
##
##tastrendh_w = np.nanmean(slopetash[0:3],axis=0)
##tastrendh_sp = np.nanmean(slopetash[3:6],axis=0)
##tastrendh_su = np.nanmean(slopetash[6:9],axis=0)
##tastrendh_f = np.nanmean(slopetash[9:12],axis=0)
##
##tastrendf_w = np.nanmean(slopetasf[0:3],axis=0)
##tastrendf_sp = np.nanmean(slopetasf[3:6],axis=0)
##tastrendf_su = np.nanmean(slopetasf[6:9],axis=0)
##tastrendf_f = np.nanmean(slopetasf[9:12],axis=0)
#
#sith_w = np.nanmean(slopesith[:,0:3],axis=1)
#sith_sp = np.nanmean(slopesith[:,3:6],axis=1)
#sith_su = np.nanmean(slopesith[:,6:9],axis=1)
#sith_f = np.nanmean(slopesith[:,9:12],axis=1)
#
#sitf_w = np.nanmean(slopesitf[:,0:3],axis=1)
#sitf_sp = np.nanmean(slopesitf[:,3:6],axis=1)
#sitf_su = np.nanmean(slopesitf[:,6:9],axis=1)
#sitf_f = np.nanmean(slopesitf[:,9:12],axis=1)
#
#tash_w = np.nanmean(slopetash[:,0:3],axis=1)
#tash_sp = np.nanmean(slopetash[:,3:6],axis=1)
#tash_su = np.nanmean(slopetash[:,6:9],axis=1)
#tash_f = np.nanmean(slopetash[:,9:12],axis=1)
#
#tasf_w = np.nanmean(slopetasf[:,0:3],axis=1)
#tasf_sp = np.nanmean(slopetasf[:,3:6],axis=1)
#tasf_su = np.nanmean(slopetasf[:,6:9],axis=1)
#tasf_f = np.nanmean(slopetasf[:,9:12],axis=1)
#
#
#sit_w = np.append(sith_w,sitf_w,axis=0)
#sit_sp = np.append(sith_sp,sitf_sp,axis=0)
#sit_su = np.append(sith_su,sitf_su,axis=0)
#sit_f = np.append(sith_f,sitf_f,axis=0)
#
#tas_w = np.append(tash_w,tasf_w,axis=0)
#tas_sp = np.append(tash_sp,tasf_sp,axis=0)
#tas_su = np.append(tash_su,tasf_su,axis=0)
#tas_f = np.append(tash_f,tasf_f,axis=0)
#
############################################################################
############################################################################
############################################################################
#### Read in individual members
#def weightThick(var,lats,types):
#    """
#    Area weights sit array 5d [ens,year,month,lat,lon] into [ens,year,month]
#    """
#    
#    if types == 'lens':
#        sityr = np.empty((var.shape[0],var.shape[1],var.shape[2]))
#        for ens in xrange(var.shape[0]):
#            for i in xrange(var.shape[1]):
#                for j in xrange(var.shape[2]):
#                    varq = var[ens,i,j,:,:]
#                    mask = np.isfinite(varq) & np.isfinite(lats)
#                    varmask = varq[mask]
#                    areamask = np.cos(np.deg2rad(lats[mask]))
#                    sityr[ens,i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
#            
#            print 'Completed: Weighting per ensemble #%s!' % ense[ens]
#    
#    elif types == 'piomas':
#        sityr = np.empty((var.shape[0],var.shape[1]))
#        for i in xrange(var.shape[0]):
#            for j in xrange(var.shape[1]):
#                varq = var[i,j,:,:]
#                mask = np.isfinite(varq) & np.isfinite(lats)
#                varmask = varq[mask]
#                areamask = np.cos(np.deg2rad(lats[mask]))
#                sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
#     
#    print '\nCompleted: Yearly weighted SIT average!' 
#    return sityr

##### Call functions     
sitmean = weightThick(sitall,lats,'lens')
tasmean = weightThick(tasall,lats,'lens')

sitall_w = np.nanmean(sitmean[:,:,0:3],axis=2)
sitall_sp = np.nanmean(sitmean[:,:,3:6],axis=2)
sitall_su = np.nanmean(sitmean[:,:,6:9],axis=2)
sitall_f = np.nanmean(sitmean[:,:,9:12],axis=2)

tasall_w = np.nanmean(tasmean[:,:,0:3],axis=2)
tasall_sp = np.nanmean(tasmean[:,:,3:6],axis=2)
tasall_su = np.nanmean(tasmean[:,:,6:9],axis=2)
tasall_f = np.nanmean(tasmean[:,:,9:12],axis=2)

sith_w = np.nanmean(sitall_w[:,:86],axis=0)
sitf_w = np.nanmean(sitall_w[:,86:],axis=0)
tash_w = np.nanmean(tasall_w[:,:86],axis=0)
tasf_w = np.nanmean(tasall_w[:,86:],axis=0)

sith_sp = np.nanmean(sitall_sp[:,:86],axis=0)
sitf_sp = np.nanmean(sitall_sp[:,86:],axis=0)
tash_sp = np.nanmean(tasall_sp[:,:86],axis=0)
tasf_sp = np.nanmean(tasall_sp[:,86:],axis=0)

sith_su = np.nanmean(sitall_su[:,:86],axis=0)
sitf_su = np.nanmean(sitall_su[:,86:],axis=0)
tash_su = np.nanmean(tasall_su[:,:86],axis=0)
tasf_su = np.nanmean(tasall_su[:,86:],axis=0)

sith_f = np.nanmean(sitall_f[:,:86],axis=0)
sitf_f = np.nanmean(sitall_f[:,86:],axis=0)
tash_f = np.nanmean(tasall_f[:,:86],axis=0)
tasf_f = np.nanmean(tasall_f[:,86:],axis=0)

### Plot PC time series
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

fig = plt.figure()
ax = plt.subplot(221)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitall_w.shape[0]):
    plt.scatter(sitall_w[i],tasall_w[i],color='darkgrey',edgecolor='darkgrey',
                s=0.05,alpha=1)
            
plt.scatter(sith_w,tash_w,color='teal',edgecolor='teal',
            s=0.7,alpha=1)
plt.scatter(sitf_w,tasf_w,color='indianred',edgecolor='indianred',
            s=0.7,alpha=1)            
#plt.plot(timex,line1,zorder=1,linewidth=2,color='indianred')
#plt.plot(smoothed1[:,0],smoothed1[:,1],color='darkslateblue',zorder=8,linewidth=1.2)

plt.xticks(np.arange(0,4,0.5),
           map(str,np.arange(0,4,0.5)),fontsize=8)
plt.xlim([0,3.5])

plt.yticks(np.arange(-34,-15,3),map(str,np.arange(-34,-15,3)),fontsize=8)
plt.ylim([-34,-16])

plt.ylabel(r'\textbf{2 m Temperature ($^\circ$C)}')

plt.text(-0.05,-16,r'\textbf{JFM}',fontsize=20,color='darkgrey')

###########################################################################

ax = plt.subplot(222)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitall_w.shape[0]):
    plt.scatter(sitall_sp[i],tasall_sp[i],color='darkgrey',edgecolor='darkgrey',
                s=0.05,alpha=1)

plt.scatter(sith_sp,tash_sp,color='teal',edgecolor='teal',
            s=0.7,alpha=1)
plt.scatter(sitf_sp,tasf_sp,color='indianred',edgecolor='indianred',
            s=0.7,alpha=1) 
#plt.plot(timex,line2,zorder=1,linewidth=2,color='indianred')
#plt.plot(smoothed2[:,0],smoothed2[:,1],color='darkslateblue',zorder=8,linewidth=1.2)

plt.xticks(np.arange(0,4,0.5),
           map(str,np.arange(0,4,0.5)),fontsize=8)
plt.xlim([0,3.5])

plt.yticks(np.arange(-15,-4,5),map(str,np.arange(-15,-4,5)),fontsize=8)
plt.ylim([-15,-5])

plt.text(2.5,-5,r'\textbf{1920-2005}',color='teal')
plt.text(2.5,-6,r'\textbf{2006-2080}',color='indianred')

plt.text(-0.05,-5,r'\textbf{AMJ}',fontsize=20,color='darkgrey')

###########################################################################

ax = plt.subplot(223)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitall_w.shape[0]):
    plt.scatter(sitall_su[i],tasall_su[i],color='darkgrey',edgecolor='darkgrey',
                s=0.05,alpha=1)

plt.scatter(sith_su,tash_su,color='teal',edgecolor='teal',
            s=0.7,alpha=1)
plt.scatter(sitf_su,tasf_su,color='indianred',edgecolor='indianred',
            s=0.7,alpha=1) 
#plt.plot(timex,line3,zorder=1,linewidth=2,color='indianred')
#plt.plot(smoothed3[:,0],smoothed3[:,1],color='darkslateblue',zorder=8,linewidth=1.2)

plt.xticks(np.arange(0,4,0.5),
           map(str,np.arange(0,4,0.5)),fontsize=8)
plt.xlim([0,3.5])

plt.yticks(np.arange(-6,7,3),map(str,np.arange(-6,7,3)),fontsize=8)
plt.ylim([-6,6])

plt.xlabel(r'\textbf{Sea Ice Thickness (m)}')
plt.ylabel(r'\textbf{2 m Temperature ($^\circ$C)}')

plt.text(-0.05,6,r'\textbf{JAS}',fontsize=20,color='darkgrey')

###########################################################################

ax = plt.subplot(224)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

for i in xrange(sitall_w.shape[0]):
    plt.scatter(sitall_f[i],tasall_f[i],color='darkgrey',edgecolor='darkgrey',
                s=0.05,alpha=1)

plt.scatter(sith_f,tash_f,color='teal',edgecolor='teal',
            s=0.7,alpha=1)
plt.scatter(sitf_f,tasf_f,color='indianred',edgecolor='indianred',
            s=0.7,alpha=1) 
#plt.plot(timex,line4,zorder=1,linewidth=2,color='indianred')
#plt.plot(smoothed4[:,0],smoothed4[:,1],color='darkslateblue',zorder=8,linewidth=1.2)

plt.xticks(np.arange(0,4,0.5),
           map(str,np.arange(0,4,0.5)),fontsize=8)
plt.xlim([0,3.5])

plt.yticks(np.arange(-30,0,5),map(str,np.arange(-30,0,5)),fontsize=8)
plt.ylim([-30,-5])

plt.text(-0.05,-5,r'\textbf{OND}',fontsize=20,color='darkgrey')

plt.xlabel(r'\textbf{Sea Ice Thickness (m)}')


fig.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'LENS_SITT2M_scatter_all.png',dpi=300)


###########################################################################
###########################################################################
###########################################################################
### Plot PC time series
### Adjust axes in time series plots 
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
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#
#fig = plt.figure()
#ax = plt.subplot(221)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendh_w,tastrendh_w,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line1,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed1[:,0],smoothed1[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,2.6,0.5),map(str,np.arange(-0.5,2.6,0.5)),fontsize=8)
#plt.ylim([-0.5,2.5])
#
#plt.text(-2.05,2.35,r'\textbf{JFM}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(222)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendh_sp,tastrendh_sp,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line2,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed2[:,0],smoothed2[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,2.6,0.5),map(str,np.arange(-0.5,2.6,0.5)),fontsize=8)
#plt.ylim([-0.5,2.5])
#
#plt.text(-2.05,2.35,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(223)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendh_su,tastrendh_su,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line3,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed3[:,0],smoothed3[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,2.6,0.5),map(str,np.arange(-0.5,2.6,0.5)),fontsize=8)
#plt.ylim([-0.5,2.5])
#
#plt.text(-2.05,2.35,r'\textbf{JAS}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(224)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendh_f,tastrendh_f,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line4,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed4[:,0],smoothed4[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,2.6,0.5),map(str,np.arange(-0.5,2.6,0.5)),fontsize=8)
#plt.ylim([-0.5,2.5])
#
#plt.text(-2.05,2.35,r'\textbf{OND}',fontsize=20,color='darkgrey')
#
#
#fig.subplots_adjust(hspace=0.4)
#plt.savefig(directoryfigure + 'LENS_SITT2Mtrends_scatter.png',dpi=300)
#
#
#
#fig = plt.figure()
#ax = plt.subplot(221)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendf_w,tastrendf_w,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line1,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed1[:,0],smoothed1[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,3.6,0.5),map(str,np.arange(-0.5,3.6,0.5)),fontsize=8)
#plt.ylim([-0.5,3.5])
#
#plt.text(-2.05,3.4,r'\textbf{JFM}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(222)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendf_sp,tastrendf_sp,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line2,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed2[:,0],smoothed2[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,3.6,0.5),map(str,np.arange(-0.5,3.6,0.5)),fontsize=8)
#plt.ylim([-0.5,3.5])
#
#plt.text(-2.05,3.4,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(223)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendf_su,tastrendf_su,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line3,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed3[:,0],smoothed3[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,3.6,0.5),map(str,np.arange(-0.5,3.6,0.5)),fontsize=8)
#plt.ylim([-0.5,3.5])
#
#plt.text(-2.05,3.4,r'\textbf{JAS}',fontsize=20,color='darkgrey')
#
############################################################################
#
#ax = plt.subplot(224)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#plt.scatter(sittrendf_f,tastrendf_f,color='k',edgecolor='k',
#            s=0.3,alpha=0.5)
##plt.plot(timex,line4,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed4[:,0],smoothed4[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(-2,2.1,0.5),
#           map(str,np.arange(-2,2.1,0.5)),fontsize=8)
#plt.xlim([-2,2])
#
#plt.yticks(np.arange(-0.5,3.6,0.5),map(str,np.arange(-0.5,3.6,0.5)),fontsize=8)
#plt.ylim([-0.5,3.5])
#
#plt.text(-2.05,3.4,r'\textbf{OND}',fontsize=20,color='darkgrey')
#
#
#fig.subplots_adjust(hspace=0.4)
#plt.savefig(directoryfigure + 'LENS_SITT2Mtrends_scatter_rcp85.png',dpi=300)
