"""
Script plots relationships between SIV and SIT (1979-2015)
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 15 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap
import operator
import matplotlib.colors as c
from calc_SIV import sivYear
import scipy.stats as sts

### Define directories
directorydata = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Calculate sea ice volume - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

### Call functions
lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
lats,lons,sic = CC.readPiomas(directorydata,years,0.01)
area = CA.readPiomasArea(directorydata)
sivyr = sivYear(sit,sic,area,False)

def weightThick(var,area):
    """
    Area weights sit array 4d [year,month,lat,lon] into [year,month]
    """
    varyr = np.empty((var.shape[0],var.shape[1]))
    for i in xrange(var.shape[0]):
        for j in xrange(var.shape[1]):
            varq = var[i,j,:,:]
            mask = np.isfinite(varq) & np.isfinite(area)
            varmask = varq[mask]
            areamask = area[mask]
            varyr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return varyr
    
sityr = weightThick(sit,area)
sicyr = weightThick(sic,area)

### Create scatter plot
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
        
### Define colormap
def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
    
def colormapSIC():
    cmap = plt.get_cmap('RdPu')
    cmaplist = [cmap(i) for i in xrange(0,cmap.N-20)]
    cms_sic = c.ListedColormap(cmaplist)
    return cms_sic

fig = plt.figure()
ax = plt.subplot(111)

#sityr = sityr[:,0]
#sivyr = sivyr[:,0]

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

timex = np.arange(0,3.1,0.1)
varx = np.ravel(sityr)
vary = np.ravel(sivyr)
mask = np.isfinite(varx) & np.isfinite(vary)

sta = np.polyfit(varx,vary,1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx[mask],
                                                          vary[mask])
m = sta[0]
b = sta[1] 
regress = m*timex + b

plt.scatter(sityr,sivyr,c=sityr,cmap=colormapSIT())
plt.plot(timex,regress)

plt.xticks(np.arange(0,3.1,1),map(str,np.arange(0,3.1,1)))
plt.xlim([0,3])
plt.ylim([0,40])
plt.yticks(np.arange(0,41,10),map(str,np.arange(0,41,10)))
plt.xlabel('Sea Ice Thickness (m)')
plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')  

plt.grid(color='k',linewidth=0.7)
plt.text(2.6,3,r'r$^{2}$=%s' % round(r1,2))

fig.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure + 'sivsit_scatter.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

timex = np.arange(0,3.1,0.1)
varx = np.ravel(sicyr)
vary = np.ravel(sivyr)
mask = np.isfinite(varx) & np.isfinite(vary)

sta = np.polyfit(varx,vary,1)
slope,intercept,r1,p_value,std_err = sts.stats.linregress(varx[mask],
                                                          vary[mask])
m = sta[0]
b = sta[1] 
regress = m*timex + b

plt.scatter(sicyr,sivyr,c=sicyr,cmap=colormapSIC())
plt.plot(timex,regress)

plt.xticks(np.arange(0,1.1,0.25),map(str,np.arange(0,1.1,0.25)))
plt.xlim([0,1])
plt.ylim([0,40])
plt.yticks(np.arange(0,41,10),map(str,np.arange(0,41,10)))
plt.xlabel('Sea Ice Concentration (fraction)')
plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')  

plt.grid(color='k',linewidth=0.7)
plt.text(0.01,0.1,r'r$^{2}$=%s' % round(r1,2))

fig.subplots_adjust(bottom=0.15)

plt.savefig(directoryfigure + 'sicsit_scatter.png',dpi=300)