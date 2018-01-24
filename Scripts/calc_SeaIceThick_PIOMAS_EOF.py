"""
Script calcules EOFs of sea ice thickness anomalies from PIOMAS
 
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
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from eofs.standard import Eof
import nclcmaps as ncm

### Define directories
directorydatap = '/surtsey/zlabe/seaice_obs/PIOMAS/Thickness/'  
directorydata = '/surtsey/zlabe/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot SIT EOF PIOMAS - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in functions
#lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
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
sit,lat2,lon2 = readPIOMAS(directorydatap,0.15)
lats = np.unique(lat2)
lons = np.unique(lon2)

### Slice above 70
#latq = np.where(lats >=70)[0]
#lats = lats[latq]
#lats = np.squeeze(lats)
#slp = slp[:,:,latq,:]

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

meansit = climo(sit,years,1981,2010)    
anomsit = anom(meansit,sit)
    
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
eofpattern,pcpattern,varpattern = calcSeasonalEOF(anomsit,years,1979,1988,
                                       np.asarray([0,1,2,3,4,5,6,7,8,9,10,11]),2,2)

### Seasons
year1 = 1979
year2 = 2005 
years = np.arange(year1,year2+1,1) 
eof_w,pc_w,var_w = calcSeasonalEOF(anomsit,years,year1,year2,
                np.asarray([0,1,2]),2,2)
eofclimo_w,pcclimo_w,varclimo_w = calcSeasonalEOF(anomsit,years,1981,2010,
                np.asarray([0,1,2]),2,2)
eof_sp,pc_sp,var_sp = calcSeasonalEOF(anomsit,years,year1,year2,
                np.asarray([3,4,5]),2,2)
eofclimo_sp,pcclimo_sp,varclimo_sp = calcSeasonalEOF(anomsit,years,1981,2010,
                np.asarray([3,4,5]),2,2)
eof_su,pc_su,var_su = calcSeasonalEOF(anomsit,years,year1,year2,
                np.asarray([6,7,8]),2,2)
eofclimo_su,pcclimo_su,varclimo_su = calcSeasonalEOF(anomsit,years,1981,2010,
                np.asarray([6,7,8]),2,2)
eof_f,pc_f,var_f = calcSeasonalEOF(anomsit,years,year1,year2,
                np.asarray([9,10,11]),2,2)
eofclimo_f,pcclimo_f,varclimo_f = calcSeasonalEOF(anomsit,years,1981,2010,
                np.asarray([9,10,11]),2,2)                

SITindex_w = (pc_w[:,0] - np.nanmean(pcclimo_w[:,0]))/np.std(pcclimo_w[:,0])    
SITndex_sp = (pc_sp[:,0] - np.nanmean(pcclimo_sp[:,0]))/np.std(pcclimo_sp[:,0]) 
SITindex_su = (pc_su[:,0] - np.nanmean(pcclimo_su[:,0]))/np.std(pcclimo_su[:,0]) 
SITindex_f = (pc_f[:,0] - np.nanmean(pcclimo_f[:,0]))/np.std(pcclimo_f[:,0])   

var1 = eofpattern[0]
var2 = eofpattern[1]

###########################################################################
###########################################################################
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='k')
plt.rc('ytick',color='k')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='w')

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

plt.savefig(directoryfigure + 'testeofs_7905.png',dpi=300)

###########################################################################
###########################################################################