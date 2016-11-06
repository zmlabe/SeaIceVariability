"""
Script looks at NCEP/NCAR reanalysis trends
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 2 November 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_NCEP as NP
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey/NCEP/'  
directoryfigure = '/home/zlabe/Desktop/SLP/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot NCEP - %s----' % titletime 

### Alott time series
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in functions
lats,lons,slp = NP.readNCEP(directorydata,years,'slp','surface')    

latq = np.where(lats >= 70)[0]   
lats = lats[latq]
slp = slp[:,:,latq,:] 

### calculate climo
def climo(var,years,yearmin,yearmax):
    """
    Calculates climatology based on given years
    """
    yr = np.where((years >= yearmin) & (years <= yearmax))[0]
    
    meanvar = np.nanmean(var[yr,:,:,:],axis=0)
    
    print 'Completed: Calculated mean climatology!'
    return meanvar

### Calculate anomalies    
meanslp = climo(slp,years,1981,2010)
anomalies = slp - meanslp

yr = np.where(years == 2012)[0][0]

springanom = np.nanmean(anomalies[:,3:6,:,:],axis=1)
summeranom = np.nanmean(anomalies[:,6:9,:,:],axis=1)

### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

barlim = np.arange(-10,12,2)
values = np.arange(-10,10.1,0.5)
cmap = ncm.cmap('testcmap')       

ax1 = plt.subplot(611)
plt.contourf(anomalies[yr,3,:,:],values,cmap=cmap,extend='both')

ax1.spines['top'].set_color('darkgrey')
ax1.spines['right'].set_color('darkgrey')
ax1.spines['bottom'].set_color('darkgrey')
ax1.spines['left'].set_color('darkgrey')
plt.setp(ax1.get_yticklabels(), visible=False)
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.xaxis.set_tick_params(size=0)
ax1.yaxis.set_tick_params(size=0)
plt.grid(color='darkgrey')
ax1.yaxis.grid(False) 

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax1.annotate(r'\textbf{A}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

ax2 = plt.subplot(612)
plt.contourf(anomalies[yr,4,:,:],values,cmap=cmap,extend='both')

ax2.spines['top'].set_color('darkgrey')
ax2.spines['right'].set_color('darkgrey')
ax2.spines['bottom'].set_color('darkgrey')
ax2.spines['left'].set_color('darkgrey')
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
ax2.xaxis.set_tick_params(size=0)
ax2.yaxis.set_tick_params(size=0)
plt.grid(color='darkgrey')
ax2.yaxis.grid(False) 

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax2.annotate(r'\textbf{M}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

ax3 = plt.subplot(613)
plt.contourf(anomalies[yr,5,:,:],values,cmap=cmap,extend='both')

ax3.spines['top'].set_color('darkgrey')
ax3.spines['right'].set_color('darkgrey')
ax3.spines['bottom'].set_color('darkgrey')
ax3.spines['left'].set_color('darkgrey')
plt.setp(ax3.get_yticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
ax3.xaxis.set_tick_params(size=0)
ax3.yaxis.set_tick_params(size=0)
plt.grid(color='darkgrey')
ax3.yaxis.grid(False) 

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax3.annotate(r'\textbf{J}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

ax4 = plt.subplot(614)
plt.contourf(anomalies[yr,6,:,:],values,cmap=cmap,extend='both')

ax4.spines['top'].set_color('darkgrey')
ax4.spines['right'].set_color('darkgrey')
ax4.spines['bottom'].set_color('darkgrey')
ax4.spines['left'].set_color('darkgrey')
plt.setp(ax4.get_yticklabels(), visible=False)
plt.setp(ax4.get_xticklabels(), visible=False)
ax4.xaxis.set_tick_params(size=0)
ax4.yaxis.set_tick_params(size=0)
plt.grid(color='darkgrey')
ax4.yaxis.grid(False) 

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax4.annotate(r'\textbf{J}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

ax5 = plt.subplot(615)
plt.contourf(anomalies[yr,7,:,:],values,cmap=cmap,extend='both')

ax5.spines['top'].set_color('darkgrey')
ax5.spines['right'].set_color('darkgrey')
ax5.spines['bottom'].set_color('darkgrey')
ax5.spines['left'].set_color('darkgrey')
plt.setp(ax5.get_yticklabels(), visible=False)
plt.setp(ax5.get_xticklabels(), visible=False)
ax5.xaxis.set_tick_params(size=0)
ax5.yaxis.set_tick_params(size=0)
plt.grid(color='darkgrey')
ax5.yaxis.grid(False) 

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax5.annotate(r'\textbf{A}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

ax6 = plt.subplot(616)
cs = plt.contourf(anomalies[yr,8,:,:],values,cmap=cmap,extend='both')

ax6.spines['top'].set_color('darkgrey')
ax6.spines['right'].set_color('darkgrey')
ax6.spines['left'].set_color('darkgrey')
ax6.spines['bottom'].set_color('darkgrey')
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off')    
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on')   
plt.grid(color='darkgrey')
ax6.yaxis.grid(False) 

ax6.annotate(r'\textbf{S}', xy=(0, 0), xytext=(-0.07, 0.3),
            xycoords='axes fraction',fontsize=22,color='darkgrey')

plt.xticks(np.arange(0,154,12),map(str,np.arange(0,361,30)))
ax6.get_xaxis().set_tick_params(direction='out', width=1,length=6,color='darkgrey')
plt.xlim([0,144])

fig.subplots_adjust(hspace=0.02)
fig.subplots_adjust(bottom=0.19)

cbar_ax = fig.add_axes([0.215,0.06,0.6,0.02])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)
#cbar.set_label(r'\textbf{H7( m )}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 
cbar.ax.tick_params(axis='x', size=.1)

plt.annotate(r'\textbf{SLP( mb )}', xy=(0, 0), xytext=(0.48,0.09),
            xycoords='figure fraction',fontsize=9,color='k')


plt.savefig(directoryfigure + 'slp_seasons_2012.png',
            dpi=300)
print 'Completed: Script done!'