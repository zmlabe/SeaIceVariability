"""
Scripts calculates Arctic Dipole using 2nd eof of MSLP anomaly north
of 70N. 
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 17 February 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_20CRv2 as CR
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from eofs.standard import Eof
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey/20CR_v2/'  
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
yearmin = 1851
yearmax = 2014
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in functions
lats,lons,slp = CR.read20CR(directorydata,years,'prmsl','surface') 

### Slice above 70
latq = np.where(lats >=70)[0]
lats = lats[latq]
lats = np.squeeze(lats)
slp = slp[:,:,latq,:]

### calculate climo
def climo(var,years,yearmin,yearmax):
    """
    Calculates climatology based on given years
    """
    print '\n>>> Using climo function!'
    yr = np.where((years >= yearmin) & (years <= yearmax))[0]
    
    meanvar = np.nanmean(var[yr,:,:,:],axis=0)
    
    print '*Completed: Calculated mean climatology!'
    return meanvar

### Calculate anomalies  
def anom(meanslp,slp):
    """
    Calculates SLP anomalies
    """
    print '\n>>> Using anom function!'
    
    anomslp = slp - meanslp
    print 'Completed: Calculated anomalies!'
    return anomslp

meanslp = climo(slp,years,1981,2010)    
anomslp = anom(meanslp,slp)
    
def calcSeasonalEOF(anomslp,years,year1,year2,monthind,eoftype,pctype):
    """
    Calculates EOF over defined seasonal period
    
    Parameters
    ----------
    anomslp : 4d array [year,month,lat,lon]
        sea level pressure anomalies
    years : 1d array
        years in total
    year1 : integer
        min month
    year2 : integer
        max month
    monthind : 1d array
        indices for months to be calculated in seasonal mean
    eoftype : integer
        1,2
    pctype : integer
        1,2
    
    Returns
    -------
    eof : array
        empirical orthogonal function
    pc : array
        principal components
    """
    print '\n>>> Using calcSeasonalEOF function!'
    
    ### Slice years
    if np.isfinite(year1):
        if np.isfinite(year2):
            yearqq = np.where((years >= year1) & (years <= year2))
            anomslp = anomslp[yearqq,:,:,:].squeeze()
        else:
            print 'Using entire time series for this EOF!'
    else:
        print 'Using entire time series for this EOF!'   
    print 'Sliced time period for seasonal mean!'
    
    ### Average over months
    anomslp = anomslp[:,monthind,:,:]
    anomslp = np.nanmean(anomslp[:,:,:,:],axis=1)
    
    print 'Sliced month period for seasonal mean!'
    
    ### Calculate EOF
    # Create an EOF solver to do the EOF analysis. Square-root of cosine of
    # latitude weights are applied before the computation of EOFs.
    coslat = np.cos(np.deg2rad(lats)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    solver = Eof(anomslp, weights=wgts)
    
    # Retrieve the leading EOF, expressed as the covariance between the 
    # leading PC time series and the input SLP anomalies at each grid point.
    eof = solver.eofsAsCovariance(neofs=eoftype)
    pc = solver.pcs(npcs=pctype, pcscaling=1)
    
    print 'EOF and PC computed!'
    
    print '*Completed: EOF and PC Calculated!\n'
    
    return eof,pc
    
### Climatology
eofpattern,pcpattern = calcSeasonalEOF(anomslp,years,1981,2010,
                                       np.asarray([0,1,2,10,11]),2,2)
                                       
if eofpattern[1,3,-10] < 0:
    eofpattern[1,:,:] *=-1                                      

### Seasons
year1 = 1851
year2 = 2014  
years2 = np.arange(year1,year2+1,1) 
eof_w,pc_w = calcSeasonalEOF(anomslp,years,year1,year2,
                np.asarray([0,1,2]),2,2)
eofclimo_w,pcclimo_w = calcSeasonalEOF(anomslp,years,1851,2014,
                np.asarray([0,1,2]),2,2)
eof_sp,pc_sp = calcSeasonalEOF(anomslp,years,year1,year2,
                np.asarray([3,4,5]),2,2)
eofclimo_sp,pcclimo_sp = calcSeasonalEOF(anomslp,years,1851,2014,
                np.asarray([3,4,5]),2,2)
eof_su,pc_su = calcSeasonalEOF(anomslp,years,year1,year2,
                np.asarray([6,7,8]),2,2)
eofclimo_su,pcclimo_su = calcSeasonalEOF(anomslp,years,1851,2014,
                np.asarray([6,7,8]),2,2)
eof_f,pc_f = calcSeasonalEOF(anomslp,years,year1,year2,
                np.asarray([9,10,11]),2,2)
eofclimo_f,pcclimo_f = calcSeasonalEOF(anomslp,years,1851,2014,
                np.asarray([9,10,11]),2,2)                

AOindex_w = (pc_w[:,0] - np.nanmean(pcclimo_w[:,0]))/np.std(pcclimo_w[:,0])    
AOindex_sp = (pc_sp[:,0] - np.nanmean(pcclimo_sp[:,0]))/np.std(pcclimo_sp[:,0]) 
AOindex_su = (pc_su[:,0] - np.nanmean(pcclimo_su[:,0]))/np.std(pcclimo_su[:,0]) 
AOindex_f = (pc_f[:,0] - np.nanmean(pcclimo_f[:,0]))/np.std(pcclimo_f[:,0])               
                                       
### Calculate Dipole
def ArcticDipoleIndex(anomslp,eofpattern,years,year1,year2,monthind):
    """
    Calculate Arctic Dipole index using Overland et al. [2010]
    convention by regressing slp anomalies to EOF spatial pattern
    and then averaging per season
    """
    print '\n>>> Using Arctic DipoleIndex function!'       

    ### Slice years
    if np.isfinite(year1):
        if np.isfinite(year2):
            yearqq = np.where((years >= year1) & (years <= year2))
            anomslp = anomslp[yearqq,:,:,:].squeeze()
        else:
            print 'Using entire time series for this EOF!'
    else:
        print 'Using entire time series for this EOF!'   
    print 'Sliced time period for seasonal mean!'  
    
    ad = np.empty((anomslp.shape[0],anomslp.shape[1]))
    for i in xrange(anomslp.shape[0]):
        print 'Regressing year ---> %s!' % (np.arange(year1,year2+1,1)[i])
        for j in xrange(anomslp.shape[1]):
            varx = np.ravel(anomslp[i,j,:,:])
            vary = np.ravel(eofpattern[1,:,:])
            mask = np.isfinite(varx) & np.isfinite(vary)     
            
            ad[i,j],intercept,r,p_value,std_err = sts.stats.linregress(varx[mask],
                                                                  vary[mask]) 
                                                                  
    ### Average over months
    ad = ad[:,monthind]
    ad = np.nanmean(ad[:,:],axis=1)     
                                                         
    print '*Completed: AD Regressed!\n'                                                             
    return ad

### Seasons
ad_w = ArcticDipoleIndex(anomslp,eofpattern,years,year1,
                       year2,np.asarray([0,1,2]))  
adclimo_w = ArcticDipoleIndex(anomslp,eofpattern,years,1851,
                       2014,np.asarray([0,1,2]))   
ad_sp = ArcticDipoleIndex(anomslp,eofpattern,years,year1,
                       year2,np.asarray([3,4,5]))  
adclimo_sp = ArcticDipoleIndex(anomslp,eofpattern,years,1851,
                       2014,np.asarray([3,4,5]))  
ad_su = ArcticDipoleIndex(anomslp,eofpattern,years,year1,
                       year2,np.asarray([6,7,8]))  
adclimo_su = ArcticDipoleIndex(anomslp,eofpattern,years,1851,
                       2014,np.asarray([6,7,8]))  
ad_f = ArcticDipoleIndex(anomslp,eofpattern,years,year1,
                       year2,np.asarray([9,10,11]))  
adclimo_f = ArcticDipoleIndex(anomslp,eofpattern,years,1851,
                       2014,np.asarray([9,10,11]))                                                             
    
adindex_w = (ad_w - np.nanmean(adclimo_w))/np.std(adclimo_w)  
adindex_sp = (ad_sp - np.nanmean(adclimo_sp))/np.std(adclimo_sp)  
adindex_su = (ad_su - np.nanmean(adclimo_su))/np.std(adclimo_su)  
adindex_f = (ad_f - np.nanmean(adclimo_f))/np.std(adclimo_f)    
                            
var1 = eofpattern[0]
var2 = eofpattern[1]

###########################################################################
###########################################################################
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 
plt.rc('axes',edgecolor='k')
plt.rc('xtick',color='k')
plt.rc('ytick',color='k')
plt.rc('axes',labelcolor='darkgrey')
plt.rc('axes',facecolor='w')

fig = plt.figure()
ax = plt.subplot(121)

m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

# Make the plot continuous
barlim = np.arange(-3,4,1)
values = np.arange(-3,3.1,0.25)

var1, lons_cyclic = addcyclic(var1, lons)
var1, lons_cyclic = shiftgrid(180., var1, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var1,
                values,extend='both')
cs1 = m.contour(x,y,var1,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-')
                
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)   

ax.annotate(r'\textbf{EOF1}',xy=(0,0),xytext=(0.35,1.05),
            textcoords='axes fraction',fontsize=20,color='darkgrey')            
                
###########################################################################                
                
ax = plt.subplot(122)
m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
            resolution='l',round =True)
m.drawmapboundary(fill_color='white')
m.drawcoastlines(color='k',linewidth=0.3)
parallels = np.arange(50,90,10)
meridians = np.arange(-180,180,30)
#m.drawparallels(parallels,labels=[False,False,False,False],
#                linewidth=0,color='k',fontsize=6)
#m.drawmeridians(meridians,labels=[True,True,False,False],
#                linewidth=0,color='k',fontsize=6)
m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')

var2, lons_cyclic = addcyclic(var2, lons)
var2, lons_cyclic = shiftgrid(180., var2, lons_cyclic, start=False)
lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
x, y = m(lon2d, lat2d)

cs = m.contourf(x,y,var2,
                values,extend='both')
cs1 = m.contour(x,y,var2,
                values,linewidths=0.2,colors='darkgrey',
                linestyles='-')                
        
cmap = ncm.cmap('nrl_sirkes')         
cs.set_cmap(cmap)      

cbar_ax = fig.add_axes([0.312,0.2,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True,cmap=cmap)

cbar.set_label(r'\textbf{hPa}',color='k')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))

ax.annotate(r'\textbf{EOF2}',xy=(0,0),xytext=(0.35,1.05),
            textcoords='axes fraction',fontsize=20,color='darkgrey')

fig.subplots_adjust(bottom=0.2)

plt.savefig(directoryfigure + 'testeofs_8110_20CR.png',dpi=300)

###########################################################################
###########################################################################
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
plt.grid(color='w',zorder=1,alpha=0.2)

zeros = [0] * len(AOindex_w)

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.plot(adindex_w,marker='o',markersize=4,zorder=2,linewidth=2,
         color='steelblue',label='PC2')
         
plt.bar(np.arange(len(AOindex_w))-0.4,AOindex_w,label='PC1',color='indianred')

plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{JFM}',fontsize=20,color='darkgrey')

ax = plt.subplot(222)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

zeros = [0] * len(AOindex_sp)

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.plot(adindex_sp,marker='o',markersize=4,zorder=2,linewidth=2,
         color='steelblue',label='PC2')
         
plt.bar(np.arange(len(AOindex_sp))-0.4,AOindex_sp,label='PC1',color='indianred')

plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{AMJ}',fontsize=20,color='darkgrey')


ax = plt.subplot(223)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

zeros = [0] * len(AOindex_su)

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.plot(adindex_su,marker='o',markersize=4,zorder=2,linewidth=2,
         color='steelblue',label='PC2')
         
plt.bar(np.arange(len(AOindex_su))-0.4,AOindex_su,label='PC1',color='indianred')

plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{JAS}',fontsize=20,color='darkgrey')


ax = plt.subplot(224)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
plt.grid(color='w',zorder=1,alpha=0.2)

zeros = [0] * len(AOindex_f)

plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)

plt.plot(adindex_f,marker='o',markersize=4,zorder=2,linewidth=2,
         color='steelblue',label='PC2')
         
plt.bar(np.arange(len(AOindex_f))-0.4,AOindex_f,label='PC1',color='indianred')

plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
         
plt.xticks(np.arange(0,year2-year1+3,5),
           map(str,np.arange(year1,year2+3,5)),fontsize=6)
plt.xlim([0,year2-year1+2])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{OND}',fontsize=20,color='darkgrey')


fig.subplots_adjust(hspace=0.4)
plt.savefig(directoryfigure + 'eof2pattern_8010_20CR.png',dpi=300)


### Create text files
directorytext = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'

np.savetxt(directorytext + 'AD_JFM_%s%s_20CR.txt' % (year1,year2),
           adindex_w)
np.savetxt(directorytext + 'AD_AMJ_%s%s_20CR.txt' % (year1,year2),
           adindex_sp)
np.savetxt(directorytext + 'AD_JAS_%s%s_20CR.txt' % (year1,year2),
           adindex_su)
np.savetxt(directorytext + 'AD_OND_%s%s_20CR.txt' % (year1,year2),
           adindex_f)  
           
np.savetxt(directorytext + 'AO_JFM_%s%s_20CR.txt' % (year1,year2),
           AOindex_w)
np.savetxt(directorytext + 'AO_AMJ_%s%s_20CR.txt' % (year1,year2),
           AOindex_sp)
np.savetxt(directorytext + 'AO_JAS_%s%s_20CR.txt' % (year1,year2),
           AOindex_su)
np.savetxt(directorytext + 'AO_OND_%s%s_20CR.txt' % (year1,year2),
           AOindex_f)            

print 'Completed: Script done!'



