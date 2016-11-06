"""
Scripts calculates SIT from LENS
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 20 October 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
import datetime
import read_SeaIceThick_LENS as LENSt
import statsmodels.api as sm
from mpl_toolkits.basemap import Basemap
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey3/'
directoryfigure = '/home/zlabe/Desktop/'
#directoryfigure = '/home/zlabe/Documents/Research/SeaIceVariability/Figures/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----LENS Mean Sea Ice Thickness - %s----' % titletime 

### Alott time series
yearmin = 400
yearmax = 2200
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

#### Call functions
#lats,lons,sit = LENSt.readLENS(directorydata,0.15)


#weights = np.cos(np.deg2rad(lats))
#sit_zonalave = np.nanmean(sit,axis=2)
#
#lons,lats = np.meshgrid(lons,lats)

def weightThick(var,lats):
    """
    Area weights sit array 4d [year,month,lat,lon] into [year,month]
    """
    sityr = np.empty((var.shape[0],var.shape[1]))
    for i in xrange(var.shape[0]):
        for j in xrange(var.shape[1]):
            varq = var[i,j,:,:]
            mask = np.isfinite(varq) & np.isfinite(lats)
            varmask = varq[mask]
            areamask = np.cos(np.deg2rad(lats[mask]))
            sityr[i,j] = np.nansum(varmask*areamask)/np.sum(areamask)
     
    print '\nCompleted: Yearly weighted SIT average!' 
    return sityr
     
sitave = weightThick(sit,lats)

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

var = sit[0,8,:,:]
var[np.where(var > 10)] = np.nan

### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

def colormapSIT():
    cmap1 = plt.get_cmap('BuPu')
    cmap2 = plt.get_cmap('RdPu_r')
    cmap3 = plt.get_cmap('gist_heat_r')
    cmaplist1 = [cmap1(i) for i in xrange(cmap1.N-10)]
    cmaplist2 = [cmap2(i) for i in xrange(15,cmap2.N)]
    cmaplist3 = [cmap3(i) for i in xrange(cmap2.N)]
    cms_sit = c.ListedColormap(cmaplist1 + cmaplist2 + cmaplist3)
    return cms_sit
    
cmap = colormapSIT()
#cmap = ncm.cmap('helix') 

#fig = plt.figure()
#ax = plt.subplot(111)
#
#m = Basemap(projection='npstere',boundinglat=62,lon_0=270,
#            resolution='l',round =True)
#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(color='k',linewidth=0.3)
#parallels = np.arange(50,90,10)
#meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0.3,color='k',fontsize=6)
#m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
#
## Make the plot continuous
#barlim = np.arange(0,8,1)
#values = np.arange(0,7.1,.25)
#
#cs = m.contourf(lons,lats,var[:,:],
#                values,extend='max',latlon=True)
#cs1 = m.contour(lons,lats,var[:,:],
#                values,linewidths=0.2,colors='k',
#                linestyles='-',latlon=True)                       
#        
#cs.set_cmap(cmap)
#
#cbar = m.colorbar(cs,location='right',pad='10%',drawedges=True)
#cbar.set_ticks(barlim)
#cbar.set_ticklabels(map(str,barlim))  
#cbar.ax.tick_params(axis='x', size=.1)
#cbar.set_label(r'\textbf{thickness (m)}')
#
#fig.suptitle(r'\textbf{LENS SIT 01/400')
#
#plt.savefig(directoryfigure + 'lens_sit_1deg.png',dpi=300)
#'Completed: Script done!'

############################################################################
############################################################################
############################################################################
sitmean = np.squeeze(np.apply_over_axes(np.nanmean,sit,[2,3]))

sitm = np.ravel(sitave[-1000:,:])
time = np.arange(0,sitm.shape[0]-1,600)

fig = plt.figure()
for i,times in enumerate(time):
    adjax = i+1  
    
    line3 = [3]*(50*12)
    
    ax = plt.subplot(5,4,i+1)
    ax.plot(line3,color='k',linewidth=0.85)
    ax.plot(np.ravel(sitm[times:times+(50*12)]),linewidth=0.7,
            color='darkslateblue')

            
    if any([adjax == 1,adjax == 5,adjax == 9,adjax == 13,adjax==17]):       
        ax.spines['top'].set_color('none')      
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)
        ax.xaxis.set_tick_params(size=0)
        ax.yaxis.set_tick_params(size=0)
        plt.setp(ax.get_xticklabels(), visible=False)
                
        plt.xlim([1,600])
        plt.ylim([2,4])
        plt.yticks(np.arange(2,5,1),map(str,np.arange(2,5,1)),fontsize=8)
        
    elif any([adjax == 4,adjax == 8,adjax == 12,adjax == 16,adjax == 20]):        
        ax.spines['top'].set_color('none')     
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        ax.xaxis.set_tick_params(size=0)
        ax.yaxis.set_tick_params(size=0)
        ax.yaxis.set_ticks_position('right')
                      
        plt.xlim([1,600])
        plt.ylim([2,4])
        plt.yticks(np.arange(2,5,1),map(str,np.arange(2,5,1)),fontsize=8)
               
    else:
        ax.spines['top'].set_color('none')     
        ax.spines['bottom'].set_linewidth(2)
        ax.spines['left'].set_linewidth(2)
        ax.spines['right'].set_linewidth(2)
        ax.xaxis.set_tick_params(size=0)
        ax.yaxis.set_tick_params(size=0)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        
        plt.xlim([1,600])
        plt.ylim([2,4])
        plt.yticks(np.arange(2,5,1),map(str,np.arange(2,5,1)),fontsize=8)
        
    plt.annotate('Period %s' % (i+1),xy=(0.5,1),xycoords='axes fraction',
                 ha='center',va='center',fontsize=6,
                 bbox=dict(facecolor='white', edgecolor='k',boxstyle='round'))
    
    plt.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'timeseries_lens.png',dpi=300)

############################################################################
############################################################################
############################################################################
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.5,width=2,which='major')  
#plt.grid(color='k',zorder=1,alpha=0.4)

plt.plot(sitave.transpose(),color='darkslateblue',alpha=0.035,linewidth=0.7,
         zorder=2)
plt.plot(np.nanmean(sitave.transpose(),axis=1),color='r',linewidth=2,
         marker='o',markeredgecolor='r',zorder=3)

plt.ylabel('Sea Ice Thickness (m)')
plt.xticks(np.arange(0,12,1),months) 
plt.yticks(np.arange(0,4.5,0.5),map(str,np.arange(0,4.5,0.5))) 
plt.xlim([0,11]) 
plt.ylim([1.5,4])

plt.savefig(directoryfigure + 'seasonalcycle_lens.png',dpi=300)

############################################################################
############################################################################
############################################################################
sitperiods = np.ravel(sitave)
sitperiod30 = np.transpose(np.reshape(sitperiods,(sitperiods.shape[0]/(30*12),(30*12))))

sitperiodmax = np.nanmax(sitperiod30,axis=1)
sitperiodmin = np.nanmin(sitperiod30,axis=1)
sitperiodmean = np.nanmean(sitperiod30,axis=1)

fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.tick_params('both',length=4.5,width=2,which='major')  
#plt.grid(color='k',zorder=1,alpha=0.4)

plt.plot(sitperiod30,color='dimgrey',alpha=0.2,
         linewidth=0.7,zorder=3)
plt.plot(sitperiodmean,color='darkslateblue',linewidth=1.5,zorder=4)
plt.plot(sitperiodmax,color='indianred',linewidth=1,zorder=2)
plt.plot(sitperiodmin,color='indianred',linewidth=1,zorder=1)

plt.xlabel(r'Years')
plt.ylabel('Sea Ice Thickness (m)')
plt.yticks(np.arange(1.5,4.5,0.5),map(str,np.arange(1.5,4.5,0.5)))         
plt.xticks(np.arange(0,361,60),map(str,np.arange(0,31,5)))
plt.xlim([0,360])
         
plt.savefig(directoryfigure + '30yearSIT_LENScontrol.png',dpi=300)