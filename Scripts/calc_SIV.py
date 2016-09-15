"""
Script calculates SIV time series 
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 12 September 2016
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

### Calculate sea ice area per grid cell    
def sivGrid(sit,sic,area,conc):
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
    sivgrid : 4d array [year,month,lat,lon]
        sea ice volume

    Usage
    -----
    sivgr = sivYear(sit,sic,area,conc)
    """
    
    print '\n>>> Using sivGrid function!'   
    
    print 'Calculating sea ice volume'
    if conc == True:
        siv = sit*sic*area
    elif conc == False:
        siv = sit*area
    else:
        RuntimeError('Did not use correct argument for function!')
        
    ### Correct units (I don't know why yet!!!!)
    siv = siv                                  
      
    print '*Completed: Calculated sea ice volume per grid cell!'
    return siv
 
###########################################################################
###########################################################################
###########################################################################   

### Call function
sivyr = sivYear(sit,sic,area,False)
#sivgr = sivGrid(sit,sic,area,False)

### Create figure for monthly SIV over time
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

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

color = iter(plt.cm.Paired(np.linspace(0,1,sivyr.shape[1])))

for i in xrange(sivyr.shape[1]):
    cma = next(color)
    plt.plot(years[:],sivyr[:,i],c=cma,linewidth=2,linestyle='-',
             label='%s' % months[i])
    l = plt.legend(shadow=False,fontsize=6,loc='upper right',
                   fancybox=True,ncol=4,bbox_to_anchor = [1.01, 1.03])
                   
    xlabels = map(str,np.arange(1979,2016,5))
    plt.xticks(np.arange(1979,2016,5),xlabels)
    plt.xlim([1979,2016])
    ylabels = map(str,np.arange(0,36,5))
    plt.yticks(np.arange(0,36,5),ylabels)
    plt.ylim([0,35])
    plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')
ax.tick_params('both',length=6,width=1.5,which='major')
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
plt.savefig(directoryfigure + 'siv_months.png',dpi=300)

### Create figure for anomalies
monthmeans = np.nanmean(sivyr,axis=0)
anoms = sivyr - monthmeans

###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
ax = plt.subplot(111)

zeroline = [0]*len(years)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)
color = iter(plt.cm.Paired(np.linspace(0,1,sivyr.shape[1])))
for i in xrange(sivyr.shape[1]):
    cma = next(color)
    plt.plot(years,anoms[:,i],c=cma,linewidth=2,linestyle='-',
             label='%s' % months[i],zorder=2)
    plt.plot(years,zeroline,linestyle='--',color='k',zorder=1,
             linewidth=2.5)
    l = plt.legend(shadow=False,fontsize=6,loc='upper right',
                   fancybox=True,ncol=4,bbox_to_anchor = [1.01, 1.03])
    xlabels = map(str,np.arange(1979,2016,5))
    plt.xticks(np.arange(1979,2016,5),xlabels)
    plt.xlim([1979,2016])
    ylabels = map(str,np.arange(-9,10,2))
    plt.yticks(np.arange(-9,10,2),ylabels)
    plt.ylim([-9,9])
    plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')

ax.tick_params('both',length=6,width=1.5,which='major')
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
plt.text(1979.2,-8.7,r'\textbf{Anomalies, 1979-2015}',fontsize=13)
plt.savefig(directoryfigure + 'siv_monthanoms',dpi=300)

###########################################################################
###########################################################################
###########################################################################

fig = plt.figure()
ax = plt.subplot(111)

vp = plt.violinplot(sivyr,showmeans=True,showmedians=False,vert=True,widths=0.6,
                    showextrema=True)

ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.tick_params(left='on',right='off',bottom='off')
ax.tick_params('both',length=6,width=1.5,which='major')
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)

for i in vp['bodies']:
    i.set_edgecolor('darkgrey')  
    i.set_facecolor('seagreen')
vp['cbars'].set_color('k')
vp['cmaxes'].set_color('k')
vp['cmins'].set_color('k')
vp['cmeans'].set_color('k')
vp['cmaxes'].set_linewidth(0.5)        
vp['cmins'].set_linewidth(0.5) 
vp['cmeans'].set_linewidth(2)
vp['cmaxes'].set_linestyle('-')        
vp['cmins'].set_linestyle('-')   

months.insert(0,'')

plt.xlim([0,13])       
plt.xticks(np.arange(0,13,1),months) 
plt.yticks(np.arange(0,36,5),map(str,(np.arange(0,36,5)))) 
plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')  

plt.savefig(directoryfigure + 'violin_siv.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
fig = plt.figure()
ax = plt.subplot(111)

adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='k',zorder=1,alpha=0.4)

color = iter(plt.cm.plasma(np.linspace(0,1,sivyr.shape[0])))
for i in xrange(sivyr.shape[0]):
    cma = next(color)
    plt.plot(np.arange(0,12,1),sivyr[i,:],c=cma,linewidth=2,linestyle='-',
         label='%s' % years[i],zorder=2)

    plt.plot(np.arange(0,12,1),np.nanmean(sivyr,axis=0),linewidth=3,
            color='k',linestyle='-',marker='o')
    
    plt.xlim([0,11])       
    plt.xticks(np.arange(0,12,1),months[1:]) 
    plt.yticks(np.arange(0,36,5),map(str,(np.arange(0,36,5)))) 
    plt.ylabel(r'Sea Ice Volume ($\times$1000 km$^{3}$)')
    
    l = plt.legend(shadow=False,fontsize=4,loc='upper right',
               fancybox=True,ncol=8,bbox_to_anchor = [1.05, 1.08])
              
             
plt.savefig(directoryfigure + 'siv_years.png',dpi=300)
    