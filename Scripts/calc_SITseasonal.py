"""
Scripts calculates SIT time series 4 seasons
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 6 October 2016
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
year1 = 1979
year2 = 2015
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

### Call functions
lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
area = CA.readPiomasArea(directorydata)

###########################################################################
###########################################################################
###########################################################################
### Calculating temporal sit
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
sit_w = np.nanmean(sitave[:,0:3],axis=1)
sit_sp = np.nanmean(sitave[:,3:6],axis=1)
sit_su = np.nanmean(sitave[:,6:9],axis=1)
sit_f = np.nanmean(sitave[:,9:12],axis=1)

sitz_w = (sit_w - np.nanmean(sit_w))/np.std(sit_w)
sitz_sp = (sit_sp - np.nanmean(sit_sp))/np.std(sit_sp)
sitz_su = (sit_su - np.nanmean(sit_su))/np.std(sit_su)
sitz_f = (sit_f - np.nanmean(sit_f))/np.std(sit_f)

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

fig = plt.figure()
ax = plt.subplot(221)
### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(sitz_w,marker='o',markersize=4,zorder=2,linewidth=2,
         color='forestgreen',label='SIT')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(27.5,2.7,r'\textbf{JFM}',fontsize=20)


ax = plt.subplot(222)
### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(sitz_sp,marker='o',markersize=4,zorder=2,linewidth=2,
         color='forestgreen',label='SIT')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(27.5,2.7,r'\textbf{AMJ}',fontsize=20)


ax = plt.subplot(223)
### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(sitz_su,marker='o',markersize=4,zorder=2,linewidth=2,
         color='forestgreen',label='SIT')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(27.5,2.7,r'\textbf{JAS}',fontsize=20)


ax = plt.subplot(224)
### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(sitz_f,marker='o',markersize=4,zorder=2,linewidth=2,
         color='forestgreen',label='SIT')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(27.5,2.7,r'\textbf{OND}',fontsize=20)

fig.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'seasonal_sit.png',dpi=300)

### Create text files
directorytext = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'

np.savetxt(directorytext + 'sit_JFM_%s%s.txt' % (year1,year2),
           sitz_w)
np.savetxt(directorytext + 'sit_AMJ_%s%s.txt' % (year1,year2),
           sitz_sp)
np.savetxt(directorytext + 'sit_JAS_%s%s.txt' % (year1,year2),
           sitz_su)
np.savetxt(directorytext + 'sit_OND_%s%s.txt' % (year1,year2),
           sitz_f)  
           
print 'Completed: Script done!'           