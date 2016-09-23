"""
Scripts creates subplots of large grid cells (nxn) for different statistical
variables. 

Author : Zachary Labe
Date : 21 September 2016
"""

### Import modules
import numpy as np
from netCDF4 import Dataset
import scipy.stats as sts
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import datetime
import iris as ir
import iris.quickplot as qplt

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/Thickness/'  
directoryfigure = '/home/zlabe/Desktop/GridCells/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate PIOMAS large grid cells - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
timex = np.arange(0,37,1)
timey = np.arange(0,8,1)

### Read in 100km EASE Piomas regridded
data = Dataset(directorydata + 'piomas_regrid_sit_19792015.nc')
lats = data.variables['lat'][:]
lons = data.variables['lon'][:]
sit = data.variables['newthickness'][:]
data.close()

sit[np.where(sit < 0.01)] = np.nan

print 'Completed: Read PIOMAS data!'

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

def transformGrid(var,la,lo):
    """
    Creates new grid with filled cells for averaged thickness
    over set bounds.
    """
    var = np.nanmean(var,axis=1)

    varn_re = np.empty(var.shape)
    for i in xrange(var.shape[0]):
        for j in xrange(0,var.shape[1]-la,la):
            for k in xrange(0,var.shape[2]-lo,lo):
                averaging = np.nanmean(var[i,j:j+la,k:k+lo])
                varn_re[i,j:j+la,k:k+lo] = averaging
                
    print 'Completed: Grid transformation!'                   
    return varn_re

la = 5
lo = 6

sitq = transformGrid(sit,la,lo)
sitq[np.where(sitq < 0.05)] = np.nan

plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

r = np.zeros((sitq.shape[1],sitq.shape[2]))
slopesit = np.zeros((sitq.shape[1],sitq.shape[2]))
intercept = np.zeros((sitq.shape[1],sitq.shape[2]))
for i in xrange(0,sitq.shape[1]-la,la):
    for j in xrange(0,sitq.shape[2]-lo,lo):
        varyy = np.ravel(sitq[:,i,j])
        varxx = np.arange(varyy.shape[0])
        mask = np.isfinite(varxx) & np.isfinite(varyy)
        
        
        if np.nanmean(sitq[:,i,j]) >= 0.15:
            fig = plt.figure()
            ax = plt.subplot(111)
            plt.scatter(timex,sitq[:,i,j],
                        color='seagreen',s=35,zorder=3,
                        edgecolor='darkgreen',linewidth=0.5)
        
            varyymean = np.nanmean(varyy)
            if np.isfinite(varyymean):
                slopesit[i:i+la,j:j+lo],intercept[i:i+la,j:j+lo],r[i:i+la,j:j+lo],p_value,std_err = sts.stats.linregress(varxx[mask],
                                                                  varyy[mask])
                fit = np.polyfit(varxx[mask],varyy[mask],1)
                m = fit[0]
                b = fit[1]
                line = m*timex + b                                
                plt.plot(line,linewidth=1,color='r',zorder=2)  
                ax.text(0.2,6.6,r'r$^2$= %s' % abs(round(r[i:i+la,j:j+lo][0][0],2)),
                        color='k',fontsize=15)  
                ax.text(0.2,6.1,r'trend = %s (cm/yr)' % (round(slopesit[i:i+la,j:j+lo][0][0],3)*100.),
                        color='k',fontsize=15)                                              
            else:
                slopesit[i:i+la,j:j+lo] = np.nan   
                r[i:i+la,j:j+lo] = np.nan
                intercept[i:i+la,j:j+lo] = np.nan
                
            plt.xlim([0,36])
            plt.ylim([0,7])
            plt.xticks(np.arange(0,37,3),map(str,np.arange(1979,2017,3)))
            plt.yticks(np.arange(0,8,1),map(str,np.arange(0,8,1)))
            plt.ylabel(r'\textbf{Thickness (m)}',fontsize=11)
            ax.spines['top'].set_color('none')
            ax.spines['right'].set_color('none')
            adjust_spines(ax, ['left', 'bottom'])
            ax.tick_params(labeltop='off', labelright='off')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            
            a2 = plt.axes([.6, .57, .29, .29], axisbg='w')   
            m = Basemap(projection='npstere',boundinglat=60,lon_0=-90,
                        resolution='l',round=True)
            m.drawmapboundary(fill_color = 'white')
            m.drawcoastlines(color = 'darkgrey',linewidth=0.2)
            m.drawlsmask(land_color='darkgrey',ocean_color='snow')
            parallels = np.arange(50,90,10)
            meridians = np.arange(-180,180,30)
            m.drawparallels(parallels,labels=[False,False,False,False],
                            linewidth=0.25)
            m.drawmeridians(meridians,labels=[True,True,True,True],
                            linewidth=0.25,
                            fontsize=4)
                            
            filled = sitq[0,i:i+la,j:j+lo].copy()
            filled.fill(1)
            cs = m.contourf(lons[i:i+la,j:j+lo],lats[i:i+la,j:j+lo],
                            filled,np.arange(0,8,1),
                            latlon=True,colors='seagreen')    
            m.fillcontinents(color='darkgrey')
            
            fig.subplots_adjust(bottom=0.15)
            
            plt.savefig(directoryfigure + '%s_%s_%s.png' % ('piomas',i,j),dpi=300)
                                      
print 'Completed: Script done!'