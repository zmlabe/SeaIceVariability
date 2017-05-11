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
#sit,lat2,lon2 = lens.readLENSEnsemble(directorydatal,0.15,'historical')
lats = np.unique(lat2)
lons = np.unique(lon2)

def readPIOMAS(directorydata,threshold):
    files = 'piomas_regrid_sit_LENS_19792015.nc'
    filename = directorydata + files
    
    data = Dataset(filename)
    sitp = data.variables['sit'][:,:,156:180,:] # lats > 65
    lats = data.variables['lat'][156:180,:]
    lons = data.variables['lon'][156:180,:]
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
    return sitp,lats,lons
    
### Read in functions
sitp,lat2,lon2 = readPIOMAS(directorydatap,0.15)
lats = np.unique(lat2)
lons = np.unique(lon2)

### calculate climo
def climo(var,years,yearmin,yearmax):
    """
    Calculates climatology based on given years
    """
    print '\n>>> Using climo function!'
    yr = np.where((years >= yearmin) & (years <= yearmax))[0]
    
    meanvar = np.nanmean(var[yr,:,:,:],axis=0)
    
    print '*Completed: Calculated mean climatology!'
    return meanvar

### Calculate anomalies  
def anom(meanslp,slp):
    """
    Calculates SLP anomalies
    """
    print '\n>>> Using anom function!'
    
    anomslp = slp - meanslp
    print 'Completed: Calculated anomalies!'
    return anomslp

meansitp = climo(sitp,years,1981,2010)    
anomsitp = anom(meansitp,sitp)

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
    # leading PC time series and the input SIT anomalies at each grid point.
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

sitn = np.reshape(sit,(39*86,12,24,360))
eofq,pcq,varii = calcSeasonalEOF(sitn,yearslens,1920,2005,
                                   np.asarray([0,1,2,3,4,5,6,7,8,9,10,11]),2,2)
                                   
eofp,pcp,varpii = calcSeasonalEOF(sitp,years,1979,2015,
                                       np.asarray([0,1,2,3,4,5,6,7,8,9,10,11]),2,2)
                                       
#eofq,pcq,varii = calcSeasonalEOF(sitn,yearslens,1920,2005,
#                                   np.asarray([3,4,5]),2,2)
#                                   
#eofp,pcp,varpii = calcSeasonalEOF(sitp,years,1979,2015,
#                                       np.asarray([3,4,5]),2,2)                                       
###########################################################################
###########################################################################
### Plot figure
var1 = eofq[0]
var2 = eofq[1]

varp1 = eofp[0]
varp2 = eofp[1]

fig = plt.figure()
ax = plt.subplot(221)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
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

perc = r'$\bf{\%}$'
ax.annotate(r'\textbf{EOF1}', xy=(0, 0), xytext=(-0.285, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')  
ax.annotate(r'\textbf{%s}%s' % (round(varii[0]*100,1),perc),xy=(0, 0),
            xytext=(-0.18, 0.8),xycoords='axes fraction',fontsize=8,
            color='k')             
                
########################################################################### 
                               
ax = plt.subplot(223)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
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

ax.annotate(r'\textbf{EOF2}', xy=(0, 0), xytext=(-0.285, -0.08),
            xycoords='axes fraction',fontsize=22,color='darkgrey')
ax.annotate(r'\textbf{%s}%s' % (round(varii[1]*100,1),perc), xy=(0, 0),
            xytext=(-0.18, 0.13),xycoords='axes fraction',fontsize=8,
            color='k')              

###########################################################################

ax = plt.subplot(222)

m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
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

varp1, lons_cyclic = addcyclic(varp1, lons)
varp1, lons_cyclic = shiftgrid(180., varp1, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,varp1,values,extend='both')
cs1 = m.contour(x,y,varp1,values,linewidths=0.2,colors='darkgrey',linestyles='-')
                
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)   

ax.annotate(r'\textbf{EOF1}', xy=(0, 0), xytext=(0.8, 0.9),
            xycoords='axes fraction',fontsize=22,color='darkgrey')   
ax.annotate(r'\textbf{%s}%s' % (round(varpii[0]*100,1),perc), xy=(0, 0),
            xytext=(0.95, 0.8),xycoords='axes fraction',fontsize=8,
            color='k')              
                
########################################################################### 
                               
ax = plt.subplot(224)
m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
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

varp2, lons_cyclic = addcyclic(varp2, lons)
varp2, lons_cyclic = shiftgrid(180., varp2, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,varp2,values,extend='both')
cs1 = m.contour(x,y,varp2,values,linewidths=0.2,colors='darkgrey',linestyles='-')                
        
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)      

ax.annotate(r'\textbf{EOF2}', xy=(0, 0), xytext=(0.8, -0.08),
            xycoords='axes fraction',fontsize=22,color='darkgrey')
ax.annotate(r'\textbf{%s}%s' % (round(varpii[1]*100,1),perc), xy=(0, 0),
            xytext=(0.95, 0.13),xycoords='axes fraction',fontsize=8,
            color='k')  

cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='Both',extendfrac=0.07,drawedges=True)

cbar.set_label(r'\textbf{m}',color='k')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
plt.setp(ax.get_xticklabels(),visible=False)

ax.annotate(r'\textbf{LENS}', xy=(0, 0), xytext=(-0.88, 1.045),
            xycoords='axes fraction',fontsize=15,color='k')
ax.annotate(r'\textbf{PIOMAS}', xy=(0, 0), xytext=(0.2, 1.045),
            xycoords='axes fraction',fontsize=15,color='k')

fig.subplots_adjust(top=0.95)
fig.subplots_adjust(bottom=0.2)
fig.subplots_adjust(wspace=-0.45)

plt.savefig(directoryfigure + 'eofs_lenspiomas.png',dpi=300)

###########################################################################
###########################################################################
########################################################################### 