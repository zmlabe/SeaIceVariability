"""
Scripts plots sit from LENS future
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 16 November 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_LENS as lens
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from netCDF4 import Dataset

### Define directories
directorydatal = '/surtsey/ypeings/'
directorydatap = '/surtsey/zlabe/seaice_obs/PIOMAS/Thickness/'  
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'

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
years2 = np.arange(2006,2080+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
ensemble = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,39,1)) + map(str,np.arange(101,106,1))

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

#### Call functions
sith,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'historical')
sitf,lats,lons = lens.readLENSEnsemble(directorydatal,0.15,'rcp85')
lons2,lats2 = np.meshgrid(lons,lats)

yearp1 = np.where((years >= 1980) & (years <= 1997))[0]
yearp2 = np.where((years >= 1998) & (years <= 2015))[0]
yearqh1 = np.where((years >= 1920) & (years <= 1962))[0]
yearqh2 = np.where((years >= 1963) & (years <= 2005))[0]
yearqf1 = np.where((years2 >= 2006) & (years2 <= 2042))[0]
yearqf2 = np.where((years2 >= 2043) & (years2 <= 2080))[0]

sitp = readPIOMAS(directorydatap,0.15)

### September 
sith_mo = np.squeeze(np.apply_over_axes(np.nanmean,sith[:,:,8,:,:],[0]))
sitf_mo = np.squeeze(np.apply_over_axes(np.nanmean,sitf[:,:,8,:,:],[0]))
sitall_mo = np.append(sith_mo,sitf_mo,axis=0)

sith1 = np.nanmean(sith_mo[yearqh1,:,:],axis=0)
sith2 = np.nanmean(sith_mo[yearqh2,:,:],axis=0)
sith3 = np.nanmean(sitall_mo[yearp1,:,:],axis=0)

sitf1 = np.nanmean(sitf_mo[yearqf1,:,:],axis=0)
sitf2 = np.nanmean(sitf_mo[yearqf2,:,:],axis=0)
sitf3 = np.nanmean(sitall_mo[yearp2,:,:],axis=0)

sitp_mo = sitp[:,8,:,:]

sitp1 = np.nanmean(sitp_mo[1:19],axis=0)
sitp2 = np.nanmean(sitp_mo[19:37],axis=0)

composites = [sith1,sith2,sith3,sitp1,sitf1,sitf2,sitf3,sitp2]

### Create subplots
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

fig = plt.figure()
    
for i in xrange(len(composites)):
    ax = plt.subplot(2,4,i+1)
    
    ### Select variable
    var = composites[i]
    
    m = Basemap(projection='npstere',boundinglat=66,lon_0=270,
                resolution='l',round =True)
                
    var, lons_cyclic = addcyclic(var, lons)
    var, lons_cyclic = shiftgrid(180., var, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
    x, y = m(lon2d, lat2d)      
      
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.2)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0.35,color='k',fontsize=1)
#    m.drawmeridians(meridians,labels=[False,False,False,False],
#                    linewidth=0.35,color='k',fontsize=1)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    ### Adjust maximum limits
    values = np.arange(0,7.1,0.25)  
    
    ### Plot filled contours    
    cs = m.contourf(x,y,var,
                    values,extend='max')
                    
    ### Set colormap                             
    cs.set_cmap('cubehelix')
    
cbar_ax = fig.add_axes([0.313,0.13,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)  
                    
cbar.set_ticks(np.arange(0,7.1,1))
cbar.set_ticklabels(map(str,np.arange(0,8,1)))    
cbar.set_label(r'\textbf{Thickness (meters)}')

plt.subplots_adjust(hspace=-0.3)
plt.subplots_adjust(wspace=0.1)
#plt.annotate(r'\textbf{Historical}', xy=(0, 0), xytext=(0.065, 0.778),
#            xycoords='figure fraction',fontsize=20,color='darkgrey',
#            rotation=90)
#plt.annotate(r'\textbf{RCP8.5}', xy=(0, 0), xytext=(0.065, 0.415),
#            xycoords='figure fraction',fontsize=20,color='darkgrey',
#            rotation=90)
            
plt.annotate(r'\textbf{LENS}', xy=(0, 0), xytext=(0.318, 0.87),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0,ha='center')
plt.annotate(r'\textbf{LENS}', xy=(0, 0), xytext=(0.595, 0.87),
            xycoords='figure fraction',fontsize=20,color='darkgrey',
            rotation=0,ha='center')
plt.annotate(r'\textbf{PIOMAS}', xy=(0, 0), xytext=(0.635, 0.847),
            xycoords='figure fraction',fontsize=7,color='darkgrey',
            rotation=0,ha='center')
plt.annotate(r'\textbf{PIOMAS}', xy=(0, 0), xytext=(0.72, 0.87),
            xycoords='figure fraction',fontsize=20,color='k',
            rotation=0)
            
plt.annotate(r'\textbf{1980-1997}', xy=(0, 0), xytext=(0.845, 0.815),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)
plt.annotate(r'\textbf{1998-2015}', xy=(0, 0), xytext=(0.845, 0.495),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)
plt.annotate(r'\textbf{1980-1997}', xy=(0, 0), xytext=(0.645, 0.815),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)                  
plt.annotate(r'\textbf{1998-2015}', xy=(0, 0), xytext=(0.645, 0.495),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=-40)

plt.annotate(r'\textbf{1920-1962}', xy=(0, 0), xytext=(0.174, 0.815),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)                               
plt.annotate(r'\textbf{1963-2005}', xy=(0, 0), xytext=(0.374, 0.815),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)
plt.annotate(r'\textbf{2006-2042}', xy=(0, 0), xytext=(0.174, 0.495),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)                     
plt.annotate(r'\textbf{2043-2080}', xy=(0, 0), xytext=(0.374, 0.495),
            xycoords='figure fraction',fontsize=7,color='k',
            rotation=0)  
                                
### Save figure
plt.savefig(directoryfigure +'sit_rcp_composites_september',dpi=500)



    
