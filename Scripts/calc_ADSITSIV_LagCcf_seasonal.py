"""
Scripts calculates Arctic Dipole using 2nd eof of MSLP anomaly north
of 70N. 
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 6 October 2016
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
year1 = 1979
year1AD = 1980
year2 = 2015
years = np.arange(year1,year2+1,1)
years2 = np.arange(year1AD,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in data
seasons = ['JFM','AMJ','JAS','OND']  

ad = np.empty((4,len(years)-1))
sit = np.empty((4,len(years)))
siv = np.empty((4,len(years)))
for i in xrange(len(seasons)):
    filename1 = 'AO_%s_%s%s.txt' % (seasons[i],year1AD,year2)
    filename2 = 'sit_%s_%s%s.txt' % (seasons[i],year1,year2)
    filename3 = 'siv_%s_%s%s.txt' % (seasons[i],year1,year2)
    
    ad[i,:] = np.genfromtxt(directorydata + filename1)
    sit[i,:] = np.genfromtxt(directorydata + filename2)
    siv[i,:] = np.genfromtxt(directorydata + filename3)
    
    print 'Reading data file %s!' % seasons[i]

### Call functions
lats,lons,sit = CT.readPiomas(directorydata2,years,0.15)
lats,lons,sic = CC.readPiomas(directorydata2,years,0.01)
area = CA.readPiomasArea(directorydata2)
sivgr = sivGrid(sit,sic,area,False)

### Seasonal SIV
siv_w = np.nanmean(sivgr[1:,0:3,:,:],axis=1)
siv_sp = np.nanmean(sivgr[1:,3:6,:,:],axis=1)
siv_su = np.nanmean(sivgr[1:,6:9,:,:],axis=1)
siv_f = np.nanmean(sivgr[1:,9:12,:,:],axis=1)

siv_season1 = np.append(siv_w,siv_sp,axis=0)
siv_season2 = np.append(siv_season1,siv_su,axis=0)
siv_season3 = np.append(siv_season2,siv_f,axis=0)

siv_season3 = np.reshape(siv_season3,
                         (siv_w.shape[0],4,siv_season1.shape[1],
                          siv_season1.shape[2]))

def lagcorr(x,y,lag,verbose=True):
    '''Compute lead-lag correlations between 2 time series.

    <x>,<y>: 1-D time series.
    <lag>: lag option, could take different forms of <lag>:
          if 0 or None, compute ordinary correlation and p-value;
          if positive integer, compute lagged correlation with lag
          upto <lag>;
          if negative integer, compute lead correlation with lead
          upto <-lag>;
          if pass in an list or tuple or array of integers, compute 
          lead/lag correlations at different leads/lags.

    Note: when talking about lead/lag, uses <y> as a reference.
    Therefore positive lag means <x> lags <y> by <lag>, computation is
    done by shifting <x> to the left hand side by <lag> with respect to
    <y>.
    Similarly negative lag means <x> leads <y> by <lag>, computation is
    done by shifting <x> to the right hand side by <lag> with respect to
    <y>.

    Return <result>: a (n*2) array, with 1st column the correlation 
    coefficients, 2nd column correpsonding p values.

    Currently only works for 1-D arrays.
    '''

    if len(x)!=len(y):
        raise('Input variables of different lengths.')

    #--------Unify types of <lag>-------------
    if np.isscalar(lag):
        if abs(lag)>=len(x):
            raise('Maximum lag equal or larger than array.')
        if lag<0:
            lag=-np.arange(abs(lag)+1)
        elif lag==0:
            lag=[0,]
        else:
            lag=np.arange(lag+1)    
    elif lag is None:
        lag=[0,]
    else:
        lag=np.asarray(lag)

    #-------Loop over lags---------------------
    result=[]
    if verbose:
        print '\n<lagcorr>: Computing lagged-correlations at lags:',lag

    for ii in lag:
        if ii<0:
            result.append(sts.pearsonr(x[:ii],y[-ii:]))
        elif ii==0:
            result.append(sts.pearsonr(x,y))
        elif ii>0:
            result.append(sts.pearsonr(x[ii:],y[:-ii]))

    result=np.asarray(result)
    
    result = result[:,0]


    return result

adq = np.ravel(np.transpose(ad))
sivgrn = siv_season3[:,:,:,:]

def deTrend(y):
    x = np.arange(y.shape[0])
    
    slopes = np.empty((y.shape[1],y.shape[2],y.shape[3]))
    intercepts = np.empty((y.shape[1],y.shape[2],y.shape[3]))
    for mo in xrange(y.shape[1]):
        for i in xrange(y.shape[2]):
            for j in xrange(y.shape[3]):
                mask = np.isfinite(y[:,mo,i,j])
                yy = y[:,mo,i,j]           
                
                if np.isfinite(np.nanmean(yy)):
                    slopes[mo,i,j], intercepts[mo,i,j], r_value, p_value, std_err = sts.linregress(x[mask],yy[mask])
                else:
                    slopes[mo,i,j] = np.nan
                    intercepts[mo,i,j] = np.nan
        print 'Regressed over month %s!' % (mo)
        
    print y.shape
    print slopes.shape
    print x.shape
    print intercepts.shape
    
    y_detrend = np.empty(y.shape)        
    for i in xrange(y.shape[0]):
        y_detrend[i,:,:,:] = y[i,:,:,:] - (slopes*x[i] + intercepts)
        print 'Detrended over year %s!' % (i)
     
    print 'Completed: Detrended SIV data!' 
    return y_detrend

sivgrn_dt = deTrend(sivgrn)

sivgrq = np.reshape(sivgrn_dt,(36*4,120,360))

laggn = np.empty((4,120,360))
for i in xrange(sivgrq.shape[1]):
    for j in xrange(sivgrq.shape[2]):
        varx = sivgrq[:,i,j]
        varx[np.where(np.isnan(varx))] = 0.
        vary = adq[:]
        laggn[:,i,j] = lagcorr(varx,vary,lag=3,verbose=True) 
       
### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

lagg3 = np.ma.array(laggn,mask=(laggn < 0.1) & (laggn > -0.1))


fig = plt.figure()
for i in xrange(laggn.shape[0]):
    ax = plt.subplot(2,2,i+1)
    
    lagg = lagg3[i]
    
    m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
                resolution='l',round =True)
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.2)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
#    m.drawparallels(parallels,labels=[False,False,False,False],
#                    linewidth=0.3,color='k',fontsize=6)
#    m.drawmeridians(meridians,labels=[False,False,False,False],
#                    linewidth=0.3,color='k',fontsize=6)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    # Make the plot continuous
    barlim = np.arange(-1,1.1,.5)
    values = np.arange(-1,1.1,0.1)
    
    cmap = plt.cm.seismic    
    
    cs = m.contourf(lons,lats,lagg,
                    values,latlon=True,cmap=cmap)
    cs1 = m.contour(lons,lats,lagg,
                    values,linewidths=0.2,colors='k',
                    linestyles='-',latlon=True)
            
#    cs.set_cmap('RdBu_r')
    cmap.set_bad('white',1)
    
    
    ax.annotate(r'Lag %s' % (i),xy=(0,0),xytext=(-0.1,0.94),
                textcoords='axes fraction',fontsize=9,color='darkgrey')
    
cbar_ax = fig.add_axes([0.312,0.1,0.4,0.03])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True,cmap=cmap)

cbar.set_label(r'\textbf{Correlation Coefficient}')
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim)) 

fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(wspace=-0.5)
fig.subplots_adjust(hspace=0.01)

plt.savefig(directoryfigure + 'testlag_ao_dt1_seasonal.png',dpi=400)