"""
Scripts calculates SIT time series for March and April 
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 9 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'   
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot Mean Sea Ice Thickness - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2016
years = np.arange(yearmin,yearmax+1,1)

### Call functions
lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
area = CA.readPiomasArea(directorydata)

###########################################################################
###########################################################################
###########################################################################
### Calculating temporal sit
def aveThick(sit):
    """
    Only calculates ave thickness for min and max time... no area weighting
    """
    minSep = np.squeeze(np.apply_over_axes(np.nanmean,sit[:,9,:,:],(1,2)))
    minMar = np.squeeze(np.apply_over_axes(np.nanmean,sit[:,3,:,:],(1,2)))
    
    print '\nCompleted: Calculated average for Sept/Mar (no weight)!'
    return minSep,minMar
    
def weightThick(var,area):
    """
    Area weights sit array 4d [year,month,lat,lon] into [year,month]
    """
    sityr = np.empty((var.shape[0],var.shape[1]))
    for i in xrange(var.shape[0]):
        for j in xrange(var.shape[1]):
            varq = var[i,j,:,:]
            mask = np.isfinite(varq) & np.isfinite(area)
            varmask = varq[mask]
            areamask = area[mask]
            sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr
     
sitave = weightThick(sit,area)
sitave = weightThick(sit,area)

minMar = sitave[:,2]
minSep = sitave[:,8]

### Plot Sep/Mar time series
fig = plt.figure()
ax = plt.subplot(111)
    
### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

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

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

### Plot figure
plt.plot(years,minMar,linestyle='-',marker='^',linewidth=2,
         markersize=4,color='steelblue',label=r'March',
         zorder=4)
plt.plot(years,minSep,linestyle='-',marker='o',linewidth=2,
         markersize=4,color='darkolivegreen',label='September',
         zorder=3)

### Add lowess smoothing         
smoothedMar = sm.nonparametric.lowess(minMar,years,it=0,frac=0.5)
smoothedSep = sm.nonparametric.lowess(minSep,years,it=0,frac=0.5)

### Plot lowess smoothing
plt.plot(smoothedMar[:,0],smoothedMar[:,1],color='r',zorder=1,
         label=r'LOWESS Mar')
plt.plot(smoothedSep[:,0],smoothedSep[:,1],color='r',zorder=2,
         label=r'LOWESS Sep')

### Adjust axes
xlabels = map(str,np.arange(1979,2017,5))
plt.xticks(np.arange(1979,2017,5),xlabels)
plt.xlim([1979,2017])

ylabels = map(str,np.arange(0,3.1,0.5))
plt.yticks(np.arange(0,3.1,0.5),ylabels)
plt.ylim([0,3])
plt.ylabel('Sea Ice Thickness (m)')

ax.tick_params('both',length=6,width=1.5,which='major')
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)

### Add legend
le = plt.legend(shadow=False,fontsize=8,loc='upper right',frameon=False,
                ncol=2)

### Save Figure
plt.savefig(directoryfigure + 'aveSIT_year.png',dpi=400)

###########################################################################
###########################################################################
###########################################################################
### Calculate anomalies
anomMar = minMar - np.nanmean(minMar)
anomSep = minSep - np.nanmean(minSep)

### Plot Sep/Mar time series anomalies
fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

### Plot anomalies
plt.plot(years,anomMar,zorder=2,linestyle='-',marker='^',linewidth=2,
         markersize=4,color='steelblue',label=r'March')
plt.plot(years,anomSep,zorder=3,linestyle='-',marker='o',linewidth=2,
         markersize=4,color='darkolivegreen',label=r'September')

### Add zero line
linezero = [0]*len(years)
plt.plot(years,linezero,linestyle='--',linewidth=2.5,color='k',zorder=1)

### Add lowess smoothing         
smoothedMaranom = sm.nonparametric.lowess(anomMar,years,it=0,frac=0.5)
smoothedSepanom = sm.nonparametric.lowess(anomSep,years,it=0,frac=0.5)

### Plot lowess smoothing
plt.plot(smoothedMaranom[:,0],smoothedMaranom[:,1],color='r',
         zorder=4,label=r'LOWESS Mar')
plt.plot(smoothedSepanom[:,0],smoothedSepanom[:,1],color='peru',
         zorder=5,label=r'LOWESS Sep')

### Adjust axes
xlabels = map(str,np.arange(1979,2017,5))
plt.xticks(np.arange(1979,2017,5),xlabels)
plt.xlim([1979,2017])

ylabels = map(str,np.arange(-1,1.1,0.5))
plt.yticks(np.arange(-1,1.1,0.5),ylabels)
plt.ylim([-1,1])
plt.ylabel('Sea Ice Thickness (m)')

ax.tick_params('both',length=6,width=1.5,which='major')
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)

### Add legend
le = plt.legend(shadow=False,fontsize=8,loc='upper right',frameon=False,
                ncol=2)

### Save figure
plt.savefig(directoryfigure + 'aveSIT_anoms_year.png',dpi=400)
    
