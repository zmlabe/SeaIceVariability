"""
Scripts calculates linear regression and LOWESS over Arctic Dipole from 
1948.
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 7 October 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
from calc_SIV import sivGrid
from mpl_toolkits.basemap import Basemap
import statsmodels.api as sm

### Define directories
directorydata = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'  
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot NCEP - %s----' % titletime 

### Alott time series
year1 = 1948
year2 = 2015
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in data
seasons = ['JFM','AMJ','JAS','OND']  

ad = np.empty((4,len(years)))
for i in xrange(len(seasons)):
    filename1 = 'AD_%s_%s%s.txt' % (seasons[i],year1,year2)
    ad[i,:] = np.genfromtxt(directorydata + filename1)

    print 'Reading data file %s!' % seasons[i]

timex = np.arange(len(ad[0,:]))
vary1 = ad[0,:] 
slope1,intercept1,r1,p_value,std_err = sts.stats.linregress(timex,vary1) 
line1 = slope1*timex + intercept1
vary2 = ad[1,:] 
slope2,intercept2,r2,p_value,std_err = sts.stats.linregress(timex,vary2) 
line2 = slope2*timex + intercept2
vary3 = ad[2,:] 
slope3,intercept3,r3,p_value,std_err = sts.stats.linregress(timex,vary3)
line3 = slope3*timex + intercept3 
vary4 = ad[3,:] 
slope4,intercept4,r4,p_value,std_err = sts.stats.linregress(timex,vary4) 
line4 = slope4*timex + intercept4

smoothed1 = sm.nonparametric.lowess(vary1,timex)
smoothed2 = sm.nonparametric.lowess(vary2,timex)
smoothed3 = sm.nonparametric.lowess(vary3,timex)
smoothed4 = sm.nonparametric.lowess(vary4,timex)

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

zeros = [0] * len(ad[0,:])

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.bar(np.arange(len(ad[0]))-0.4,ad[0],color='steelblue')
plt.plot(timex,line1,zorder=1,linewidth=1.3,color='r')
plt.plot(smoothed1[:,0],smoothed1[:,1],color='lime',zorder=8)

plt.axvline(np.where(years == 2007)[0],color='k',linestyle='--',
            linewidth=0.5)
plt.axvline(np.where(years == 2012)[0],color='k',linestyle='--',
            linewidth=0.5)

plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=4.5)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=4.5)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{JFM}',fontsize=20)

###########################################################################

ax = plt.subplot(222)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.bar(np.arange(len(ad[1]))-0.4,ad[1],color='steelblue')
plt.plot(timex,line2,zorder=1,linewidth=1.3,color='r')
plt.plot(smoothed2[:,0],smoothed2[:,1],color='lime',zorder=8)

plt.axvline(np.where(years == 2007)[0],color='k',linestyle='--',
            linewidth=0.5)
plt.axvline(np.where(years == 2012)[0],color='k',linestyle='--',
            linewidth=0.5)

plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=4.5)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=4.5)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{AMJ}',fontsize=20)

###########################################################################

ax = plt.subplot(223)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.bar(np.arange(len(ad[2]))-0.4,ad[2],color='steelblue')
plt.plot(timex,line3,zorder=1,linewidth=1.3,color='r')
plt.plot(smoothed3[:,0],smoothed3[:,1],color='lime',zorder=8)

plt.axvline(np.where(years == 2007)[0],color='k',linestyle='--',
            linewidth=0.5)
plt.axvline(np.where(years == 2012)[0],color='k',linestyle='--',
            linewidth=0.5)

plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=4.5)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=4.5)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{JAS}',fontsize=20)

###########################################################################

ax = plt.subplot(224)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.bar(np.arange(len(ad[3]))-0.4,ad[3],color='steelblue')
plt.plot(timex,line4,zorder=1,linewidth=1.3,color='r')
plt.plot(smoothed4[:,0],smoothed4[:,1],color='lime',zorder=8)

plt.axvline(np.where(years == 2007)[0],color='k',linestyle='--',
            linewidth=0.5)
plt.axvline(np.where(years == 2012)[0],color='k',linestyle='--',
            linewidth=0.5)

plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=4.5)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=4.5)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{OND}',fontsize=20)


fig.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'AD_timeseries_4815.png',dpi=300)