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
sith,lats,lons = lens.readLENSEnsemble(directorydataSIT,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydataSIT,0.15,'rcp85')

### Combine SIT periods
sitall = np.append(sith,sitf,axis=1)

### Read lhtfx
data = Dataset(directorydata2 + 'lens_regrid_LHFLX_19202080.nc')
tasall = data.variables['lhflx'][:]
data.close()
#          
##### 2D lat/lon arrays          
#lons,lats = np.meshgrid(lons,lats)

###########################################################################
###########################################################################
###########################################################################
### Read in individual members
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
            
            print 'Completed: Weighting per ensemble #%s!' % ense[ens]
    
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

##### Call functions     
#sitmean = weightThick(sitall,lats,'lens')
#tasmean = weightThick(tasall,lats,'lens')
#
#sitall_w = np.nanmean(sitmean[:,:,0:3],axis=2)
#sitall_sp = np.nanmean(sitmean[:,:,3:6],axis=2)
#sitall_su = np.nanmean(sitmean[:,:,6:9],axis=2)
#sitall_f = np.nanmean(sitmean[:,:,9:12],axis=2)
#
#tasall_w = np.nanmean(tasmean[:,:,0:3],axis=2)
#tasall_sp = np.nanmean(tasmean[:,:,3:6],axis=2)
#tasall_su = np.nanmean(tasmean[:,:,6:9],axis=2)
#tasall_f = np.nanmean(tasmean[:,:,9:12],axis=2)
#
#sith_w = np.nanmean(sitall_w[:,:86],axis=0)
#sitf_w = np.nanmean(sitall_w[:,86:],axis=0)
#tash_w = np.nanmean(tasall_w[:,:86],axis=0)
#tasf_w = np.nanmean(tasall_w[:,86:],axis=0)
#
#sith_sp = np.nanmean(sitall_sp[:,:86],axis=0)
#sitf_sp = np.nanmean(sitall_sp[:,86:],axis=0)
#tash_sp = np.nanmean(tasall_sp[:,:86],axis=0)
#tasf_sp = np.nanmean(tasall_sp[:,86:],axis=0)
#
#sith_su = np.nanmean(sitall_su[:,:86],axis=0)
#sitf_su = np.nanmean(sitall_su[:,86:],axis=0)
#tash_su = np.nanmean(tasall_su[:,:86],axis=0)
#tasf_su = np.nanmean(tasall_su[:,86:],axis=0)
#
#sith_f = np.nanmean(sitall_f[:,:86],axis=0)
#sitf_f = np.nanmean(sitall_f[:,86:],axis=0)
#tash_f = np.nanmean(tasall_f[:,:86],axis=0)
#tasf_f = np.nanmean(tasall_f[:,86:],axis=0)




sitw = np.nanmean(tasall[:,:,0:3],axis=2)
sitsp = np.nanmean(tasall[:,:,3:6],axis=2)
sitsu = np.nanmean(tasall[:,:,6:9],axis=2)
sitf = np.nanmean(tasall[:,:,9:12],axis=2)

sitw_w = np.nanmean(sitw,axis=0)
sitsp_sp = np.nanmean(sitsp,axis=0)
sitsu_su = np.nanmean(sitsu,axis=0)
sitf_f = np.nanmean(sitf,axis=0)



#### Plot PC time series
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
#for i in xrange(sitall_w.shape[0]):
#    plt.scatter(sitall_w[i],tasall_w[i],color='darkgrey',edgecolor='darkgrey',
#                s=0.05,alpha=1)
#            
#plt.scatter(sith_w,tash_w,color='teal',edgecolor='teal',
#            s=0.7,alpha=1)
#plt.scatter(sitf_w,tasf_w,color='indianred',edgecolor='indianred',
#            s=0.7,alpha=1)            
##plt.plot(timex,line1,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed1[:,0],smoothed1[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(0,4,0.5),
#           map(str,np.arange(0,4,0.5)),fontsize=8)
#plt.xlim([0,3.5])
#
#plt.yticks(np.arange(0,21,5),map(str,np.arange(0,21,5)),fontsize=8)
#plt.ylim([0,20])
#
#plt.ylabel(r'\textbf{LHFLX (W/m$^2$)}')
#
#plt.text(-0.05,20,r'\textbf{JFM}',fontsize=20,color='darkgrey')
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
#for i in xrange(sitall_w.shape[0]):
#    plt.scatter(sitall_sp[i],tasall_sp[i],color='darkgrey',edgecolor='darkgrey',
#                s=0.05,alpha=1)
#
#plt.scatter(sith_sp,tash_sp,color='teal',edgecolor='teal',
#            s=0.7,alpha=1)
#plt.scatter(sitf_sp,tasf_sp,color='indianred',edgecolor='indianred',
#            s=0.7,alpha=1) 
##plt.plot(timex,line2,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed2[:,0],smoothed2[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(0,4,0.5),
#           map(str,np.arange(0,4,0.5)),fontsize=8)
#plt.xlim([0,3.5])
#
#plt.yticks(np.arange(3,11,1),map(str,np.arange(3,11,1)),fontsize=8)
#plt.ylim([3,10])
#
#plt.text(2.5,9.5,r'\textbf{1920-2005}',color='teal')
#plt.text(2.5,8.85,r'\textbf{2006-2080}',color='indianred')
#
#plt.text(-0.05,10,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
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
#for i in xrange(sitall_w.shape[0]):
#    plt.scatter(sitall_su[i],tasall_su[i],color='darkgrey',edgecolor='darkgrey',
#                s=0.05,alpha=1)
#
#plt.scatter(sith_su,tash_su,color='teal',edgecolor='teal',
#            s=0.7,alpha=1)
#plt.scatter(sitf_su,tasf_su,color='indianred',edgecolor='indianred',
#            s=0.7,alpha=1) 
##plt.plot(timex,line3,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed3[:,0],smoothed3[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(0,4,0.5),
#           map(str,np.arange(0,4,0.5)),fontsize=8)
#plt.xlim([0,3.5])
#
#plt.yticks(np.arange(4,15,2),map(str,np.arange(4,15,2)),fontsize=8)
#plt.ylim([4,14])
#
#plt.xlabel(r'\textbf{Sea Ice Thickness (m)}')
#plt.ylabel(r'\textbf{LHFLX (W/m$^2$)}')
#
#plt.text(-0.05,14,r'\textbf{JAS}',fontsize=20,color='darkgrey')
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
#for i in xrange(sitall_w.shape[0]):
#    plt.scatter(sitall_f[i],tasall_f[i],color='darkgrey',edgecolor='darkgrey',
#                s=0.05,alpha=1)
#
#plt.scatter(sith_f,tash_f,color='teal',edgecolor='teal',
#            s=0.7,alpha=1)
#plt.scatter(sitf_f,tasf_f,color='indianred',edgecolor='indianred',
#            s=0.7,alpha=1) 
##plt.plot(timex,line4,zorder=1,linewidth=2,color='indianred')
##plt.plot(smoothed4[:,0],smoothed4[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
#
#plt.xticks(np.arange(0,4,0.5),
#           map(str,np.arange(0,4,0.5)),fontsize=8)
#plt.xlim([0,3.5])
#
#plt.yticks(np.arange(0,31,5),map(str,np.arange(0,31,5)),fontsize=8)
#plt.ylim([0,30])
#
#plt.text(-0.05,30,r'\textbf{OND}',fontsize=20,color='darkgrey')
#
#plt.xlabel(r'\textbf{Sea Ice Thickness (m)}')
#
#
#fig.subplots_adjust(hspace=0.4)
#plt.savefig(directoryfigure + 'LENS_SITLHTFX_scatter_all.png',dpi=300)


###########################################################################
###########################################################################
###########################################################################
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(sitw_w[-1] - sitw_w[0], lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-50,51,25)
values = np.arange(-50,51,2)

cs = m.contourf(x,y,var,values,
                extend='both')
#cs1 = m.contour(x,y,var,values,
#                linewidths=0.2,colors='k',
#                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JFM}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)

var, lons_cyclic = addcyclic(sitsp_sp[-1] - sitsp_sp[0], lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,values,extend='both'
                )
#cs1 = m.contour(x,y,var,values,
#                linewidths=0.2,colors='k',
#                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{AMJ}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(223)

m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(sitsu_su[-1] - sitsu_su[0], lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
                    
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,values,extend='both'
                )
#cs1 = m.contour(x,y,var,values,
#               linewidths=0.2,colors='k',
#                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{JAS}', xy=(0, 0), xytext=(-0.23, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

###########################################################################
###########################################################################

ax = plt.subplot(224)

m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)
            
var, lons_cyclic = addcyclic(sitf_f[-1] - sitf_f[0], lons)
var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)            
            
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.2)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

cs = m.contourf(x,y,var,values,
                extend='both')
#cs1 = m.contour(x,y,var,values,
#                linewidths=0.2,colors='k',
#                linestyles='-')
        
cmap = ncm.cmap('BlueWhiteOrangeRed')         
cs.set_cmap(cmap)
ax.annotate(r'\textbf{OND}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='Both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{Difference LHFX(W/m$^2$)}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'diff_lhfx.png',dpi=300)