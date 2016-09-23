"""
Script plots ranking on meshgrid plot for sea ice volume
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 23 September 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
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

### Calculate siv per year
def sivYear(sit,sic,area,conc):
    """
    Function calculates time series of sea ice volume per YEAR

    Parameters
    ----------
    sit : 4d array [year,month,lat,lon]
        sea ice thickness (m)
    sic : 4d array [year,month,lat,lon]
        sea ice concentration (fraction, 0-1)
    area : 2d array [lat,lon]
        area of grid cell in PIOMAS (km^2) 
    conc : boolean
        True or False ----> turning on/off using sea ice concentration

    Returns
    -------
    sivyr : 2d array [year,month]
        sea ice volume

    Usage
    -----
    siv_yr = sivYear(sit,sic,area,conc)
    """
    
    print '\n>>> Using sivYear function!'
    
    print 'Calculating sea ice volume'
    if conc == True:
        siv = sit*sic*area
    elif conc == False:
        siv = sit*area
    else:
        RuntimeError('Did not use correct argument for function!')
        
    ### Take temporal average of year
    siv = np.squeeze(np.apply_over_axes(np.nansum,
                                        siv[:,:,:,:],(2,3)))
                                        
    ### Correct units (I don't know why yet!!!!)
    siv = siv/10**6                                  
      
    print '*Completed: Calculated sea ice volume per year!'
    return siv
    
###########################################################################
###########################################################################
###########################################################################   
### Call function
sivyr = sivYear(sit,sic,area,False)
sivyr[np.where(sivyr==0.0)]=np.nan

sivyr = np.flipud(sivyr.transpose())

### Try ranking
rank = np.empty(sivyr.shape)
for i in xrange(sivyr.shape[0]):
    rank[i,:] = sts.rankdata(sivyr[i,:],method='min')

### Call parameters
plt.rcParams['text.usetex']=True
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Avant Garde'

### Plot first meshgrid
fig = plt.figure()
ax = plt.subplot(111)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.tick_params('both',length=0,width=0,which='major')

cs = plt.pcolormesh(rank,shading='faceted',edgecolor='w',
                    linewidth=0.3,vmin=1,vmax=38)

cs.set_cmap('viridis')

cbar = plt.colorbar(cs,orientation='horizontal')
cbar.set_ticks([])
cbar.ax.invert_xaxis()
cbar.set_label(r'\textbf{Sea Ice Volume Rank by Year}')

ylabels = ['D','N','O','S','A','J','J','M','A','M','F','J']
plt.yticks(np.arange(0.5,12.5,1),ylabels,ha='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=7)
plt.xticks(np.arange(0.5,37.5,3),map(str,np.arange(1979,2016,3)))
plt.xlim([0,37])

plt.text(-3,-5.3,r'Highest rank')
plt.text(34,-5.3,r'Lowest rank')

for i in xrange(rank.shape[0]):
    for j in xrange(rank.shape[1]):
        plt.text(j+0.5,i+0.5,'%s' % int(rank[i,j]),fontsize=6,
                 color='w',va='center',ha='center')

### Save figure
plt.savefig(directoryfigure + 'siv_ranks.png',dpi=300)

###########################################################################
###########################################################################
########################################################################### 
### Calculate climo
climyr = np.where((years>=1981) & (years<2010))[0]
sivave = np.nanmean(sivyr[:,climyr],axis=1)

anoms = sivyr.transpose() - sivave

### Plot second meshgrid
fig = plt.figure()
ax = plt.subplot(111)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.tick_params('both',length=0,width=0,which='major')

cs = plt.pcolormesh(np.flipud(anoms.transpose()),shading='faceted',
                    edgecolor='w',linewidth=0.3,vmin=-10,vmax=10)

cs.set_cmap('seismic_r')

cbar = plt.colorbar(cs,orientation='horizontal')
cbar.set_ticks(np.arange(-10,11,2))
cbar.set_ticklabels(map(str,np.arange(-10,11,2)))
cbar.set_label(r'\textbf{Sea Ice Volume Anomalies ($\times$1000\ km${^3}$})')

ylabels = ['D','N','O','S','A','J','J','M','A','M','F','J']
plt.yticks(np.arange(0.5,12.5,1),ylabels,ha='center')
yax = ax.get_yaxis()
yax.set_tick_params(pad=7)
plt.xticks(np.arange(0.5,37.5,3),map(str,np.arange(1979,2016,3)))
plt.xlim([0,37])

### Save figure
plt.savefig(directoryfigure + 'siv_anoms.png',dpi=300)