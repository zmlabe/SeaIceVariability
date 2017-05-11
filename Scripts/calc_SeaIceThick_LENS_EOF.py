"""
Script calcules EOFs of sea ice thickness anomalies from LENS
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Lindsay and Zhang [2005]
    Author : Zachary Labe
    Date   : 1 May 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from netCDF4 import Dataset
import read_SeaIceThick_LENS as lens
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from eofs.standard import Eof
import nclcmaps as ncm

### Define directories
directorydatap = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/' 
directorydatal = '/home/zlabe/Surtsey3/' 
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/eofs_sit/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot SIT EOF LENS - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
yearslens = np.arange(1920,2080+1,1)          
yearsclimo = np.arange(1981,2010+1,1)
ense = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))
          
### Read in functions
sit,lat2,lon2 = lens.readLENSEnsemble(directorydatal,0.15,'historical')
lats = np.unique(lat2)
lons = np.unique(lon2)

def calcSeasonalEOF(anomsit,years,year1,year2,monthind,eoftype,pctype):
    """
    Calculates EOF over defined seasonal period
    
    Parameters
    ----------
    anomsit : 4d array [year,month,lat,lon]
        sea ice thickness anomalies
    years : 1d array
        years in total
    year1 : integer
        min month
    year2 : integer
        max month
    monthind : 1d array
        indices for months to be calculated in seasonal mean
    eoftype : integer
        1,2
    pctype : integer
        1,2
    
    Returns
    -------
    eof : array
        empirical orthogonal function
    pc : array
        principal components
    """
    print '\n>>> Using calcSeasonalEOF function!'
    
    ### Slice years
    if np.isfinite(year1):
        if np.isfinite(year2):
            yearqq = np.where((years >= year1) & (years <= year2))
            anomsit = anomsit[yearqq,:,:,:].squeeze()
        else:
            print 'Using entire time series for this EOF!'
    else:
        print 'Using entire time series for this EOF!'   
    print 'Sliced time period for seasonal mean!'
    
    ### Average over months
    anomsit = anomsit[:,monthind,:,:]
    anomsit = np.nanmean(anomsit[:,:,:,:],axis=1)
    
    print 'Sliced month period for seasonal mean!'
    
    ### Calculate EOF
    # Create an EOF solver to do the EOF analysis. Square-root of cosine of
    # latitude weights are applied before the computation of EOFs.
    coslat = np.cos(np.deg2rad(lats)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(anomsit, weights=wgts)
    
    # Retrieve the leading EOF, expressed as the covariance between the 
    # leading PC time series and the input SLP anomalies at each grid point.
    eof = solver.eofsAsCovariance(neofs=eoftype)
    pc = solver.pcs(npcs=pctype, pcscaling=1)
    var = solver.varianceFraction(neigs=eoftype)
    
    print 'EOF and PC computed!'
    
    print '*Completed: EOF and PC Calculated!\n'
    
    return eof,pc,var
    
### Climatology
#eofn = []
#pcn = []    
#for i in xrange(sit.shape[0]):
#    eofq,pcq = calcSeasonalEOF(sit[i],yearslens,1920,2005,
#                                       np.asarray([0,1,2,10,11]),2,2)
#    eofn.append(eofq)
#    pcn.append(pcq)
#    
#eoff = np.asarray(eofn)
#pcc = np.asarray(pcn)

sitn = np.nanmean(sit,axis=0)
#sitn = np.reshape(sit,(39*86,12,24,360))
eofq,pcq,var = calcSeasonalEOF(sitn,yearslens,1920,2005,
                                   np.asarray([0,1,2,3,4,5,6,7,8,9,10,11]),2,2)



### Seasons
#year1 = 1920
#year2 = 2005
#eof_wn = []
#eof_spn = []
#eof_sun = []
#eof_fn = []
#eofclimo_wn = []
#eofclimo_spn = []
#eofclimo_sun = []
#eofclimo_fn = []
#pc_wn = []
#pc_spn = []
#pc_sun = []
#pc_fn = []
#pcclimo_wn = []
#pcclimo_spn = []
#pcclimo_sun = []
#pcclimo_fn = []
#for i in xrange(sit.shape[0]):
#    eof_wq,pc_wq = calcSeasonalEOF(sit[i],yearslens,year1,year2,
#                    np.asarray([0,1,2]),2,2)
#    eofclimo_wq,pcclimo_wq = calcSeasonalEOF(sit[i],yearslens,1920,2005,
#                    np.asarray([0,1,2]),2,2)
#    eof_spq,pc_spq = calcSeasonalEOF(sit[i],yearslens,year1,year2,
#                    np.asarray([3,4,5]),2,2)
#    eofclimo_spq,pcclimo_spq = calcSeasonalEOF(sit[i],yearslens,1920,2005,
#                    np.asarray([3,4,5]),2,2)
#    eof_suq,pc_suq = calcSeasonalEOF(sit[i],yearslens,year1,year2,
#                    np.asarray([6,7,8]),2,2)
#    eofclimo_suq,pcclimo_suq = calcSeasonalEOF(sit[i],yearslens,1920,2005,
#                    np.asarray([6,7,8]),2,2)
#    eof_fq,pc_fq = calcSeasonalEOF(sit[i],yearslens,year1,year2,
#                    np.asarray([9,10,11]),2,2)
#    eofclimo_fq,pcclimo_fq = calcSeasonalEOF(sit[i],yearslens,1920,2005,
#                    np.asarray([9,10,11]),2,2)  
#                    
#    if eof_spq[0,7,187] > 0:
#        eof_spq[0,:,:] *= -1
##        eof_spq = eof_spq * -1    
#    if eof_spq[1,-6,23] > 0:
#        eof_spq[1,:,:] *= -1            
#                    
#    eof_wn.append(eof_wq)
#    eof_spn.append(eof_spq)
#    eof_sun.append(eof_suq)
#    eof_fn.append(eof_fq)
#    eofclimo_wn.append(eofclimo_wq)
#    eofclimo_spn.append(eofclimo_spq)
#    eofclimo_sun.append(eofclimo_suq)
#    eofclimo_fn.append(eofclimo_fq) 
#
#    pc_wn.append(pc_wq)
#    pc_spn.append(pc_spq)
#    pc_sun.append(pc_suq)
#    pc_fn.append(pc_fq)
#    pcclimo_wn.append(pcclimo_wq)
#    pcclimo_spn.append(pcclimo_spq)
#    pcclimo_sun.append(pcclimo_suq)
#    pcclimo_fn.append(pcclimo_fq) 
#
#eof_w = np.asarray(eof_wn)
#eof_sp = np.asarray(eof_spn)
#eof_su = np.asarray(eof_sun)
#eof_f = np.asarray(eof_fn)
#eofclimo_w = np.asarray(eofclimo_wn)
#eofclimo_sp = np.asarray(eofclimo_spn)
#eofclimo_su = np.asarray(eofclimo_sun)
#eofclimo_f = np.asarray(eofclimo_fn)
#
#pc_w = np.asarray(pc_wn)
#pc_sp = np.asarray(pc_spn)
#pc_su = np.asarray(pc_sun)
#pc_f = np.asarray(pc_fn)
#pcclimo_w = np.asarray(pcclimo_wn)
#pcclimo_sp = np.asarray(pcclimo_spn)
#pcclimo_su = np.asarray(pcclimo_sun)
#pcclimo_f = np.asarray(pcclimo_fn)     
 
########################################################################### 
###########################################################################            
### Plot figure
#plt.rc('text',usetex=True)
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
#
#fig = plt.figure()
#
#for i in xrange(eoff.shape[0]):
#    
#    varf = eof_sp[i,1,:,:]    
#    
#    ax = plt.subplot(7,6,i+1)
#    
#    m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
#                resolution='l',round =True)
#    m.drawmapboundary(fill_color='white')
#    m.drawcoastlines(color='k',linewidth=0.3)
#    parallels = np.arange(50,90,10)
#    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0,color='k',fontsize=6)
#    m.drawmeridians(meridians,labels=[False,False,False,False],
#                    linewidth=0,color='k',fontsize=6)
#    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#    
#    # Make the plot continuous
#    barlim = np.arange(-1,1.1,0.5)
#    values = np.arange(-1,1.1,0.1)
#    
#    varf, lons_cyclic = addcyclic(varf, lons)
#    varf, lons_cyclic = shiftgrid(180., varf, lons_cyclic, start=False)
#    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
#    x, y = m(lon2d, lat2d)
#    
#    cs = m.contourf(x,y,varf,
#                    values,extend='both')
#    cs1 = m.contour(x,y,varf,
#                    values,linewidths=0.1,colors='darkgrey',
#                    linestyles='-')
#                    
#    ax.annotate(r'\textbf{%s}' % ense[i], xy=(0, 0), xytext=(-0.3, 0.9),
#            xycoords='axes fraction',fontsize=8,color='darkgrey',
#            rotation=0)
#            
#    cmap = ncm.cmap('nrl_sirkes')         
#    cs.set_cmap(cmap) 
#
#cbar_ax = fig.add_axes([0.51,0.16,0.31,0.027])                
#cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
#                    extend='both',extendfrac=0.07,drawedges=True)                      
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(map(str,barlim))  
#cbar.ax.tick_params(labelsize=8)   
#cbar.set_label(r'\textbf{m}')
#
#plt.subplots_adjust(wspace=-0.6)
#
#fig.suptitle(r'\textbf{EOF2, AMJ 1979-2005 (SIT) -- LENS}')
#
#plt.savefig(directoryfigure + 'eof2_LENS_2005.png',dpi=300)

###########################################################################
###########################################################################
### Plot figure
#var1 = np.nanmean(eof_sp[:,0,:,:],axis=0)
#var2 = np.nanmean(eof_sp[:,1,:,:],axis=0)
var1 = eofq[0]
var2 = eofq[1]

fig = plt.figure()
ax = plt.subplot(121)

m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

## Make the plot continuous
barlim = np.arange(-1,1.1,0.5)
values = np.arange(-1,1.1,0.1)

var1, lons_cyclic = addcyclic(var1, lons)
var1, lons_cyclic = shiftgrid(180., var1, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var1,values,extend='both')
cs1 = m.contour(x,y,var1,values,linewidths=0.2,colors='darkgrey',linestyles='-')
                
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)   

ax.annotate(r'\textbf{EOF1}',xy=(0,0),xytext=(0.35,1.05),
            textcoords='axes fraction',fontsize=20,color='darkgrey')            
                
########################################################################### 
                               
ax = plt.subplot(122)
m = Basemap(projection='npstere',boundinglat=67,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

var2, lons_cyclic = addcyclic(var2, lons)
var2, lons_cyclic = shiftgrid(180., var2, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var2,values,extend='both')
cs1 = m.contour(x,y,var2,values,linewidths=0.2,colors='darkgrey',linestyles='-')                
        
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)      

cbar_ax = fig.add_axes([0.312,0.2,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True,cmap=cmap)

cbar.set_label(r'\textbf{m}',color='k')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))

ax.annotate(r'\textbf{EOF2}',xy=(0,0),xytext=(0.35,1.05),
            textcoords='axes fraction',fontsize=20,color='darkgrey')

fig.subplots_adjust(bottom=0.2)

plt.savefig(directoryfigure + 'eofs_lensmean_2005.png',dpi=300)

###########################################################################
###########################################################################