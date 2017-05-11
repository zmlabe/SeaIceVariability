"""
Scripts calculates Arctic Dipole using 2nd eof of MSLP anomaly north
of 70N for LENS
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 28 November 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
import read_var_LENS as LV
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
from eofs.standard import Eof
import nclcmaps as ncm

### Define directories
directorydata = '/home/zlabe/Surtsey/LENS/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot LENS AD - %s----' % titletime 

ense = ['02','03','04','05','06','07','08','09'] + \
    map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))

### Alott time series
year1 = 1980
year2 = 2015
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']

yearslens = np.arange(1920,2080+1,1)          
yearsclimo = np.arange(1981,2010+1,1)
          
### Read in functions
slp,lats,lons = LV.readLENSEnsemble(directorydata,'SLP') 

### calculate climo
def climo(var,years,yearmin,yearmax):
    """
    Calculates climatology based on given years
    """
    print '\n>>> Using climo function!'
    yr = np.where((years >= yearmin) & (years <= yearmax))[0]
    
    meanvar = np.empty((var.shape[0],var.shape[2],var.shape[3],var.shape[4]))
    for i in xrange(var.shape[0]):
        meanvar[i,:,:,:] = np.nanmean(var[i,yr,:,:,:],axis=0)
    
    print '*Completed: Calculated mean climatology!'
    return meanvar

### Calculate anomalies  
def anom(meanslp,slp):
    """
    Calculates SLP anomalies
    """
    print '\n>>> Using anom function!'
    
    anomslp = np.empty(slp.shape)
    for i in xrange(slp.shape[0]):
        anomslp[i,:,:,:,:] = slp[i,:,:,:,:] - meanslp[i,:,:,:]
    print 'Completed: Calculated anomalies!'
    return anomslp
    
meanslp = climo(slp,yearslens,1981,2010)    
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

###########################################################################
###########################################################################
###########################################################################

eofn = []
pcn = []    
for i in xrange(slp.shape[0]):
    eofq,pcq = calcSeasonalEOF(anomslp[i],yearslens,1981,2010,
                                       np.asarray([0,1,2,10,11]),2,2)
    
#    if eofq[1,3,-10] < 0:
#        eofq = eofq * -1
    
    if eofq[0,3,-10] > 0:
        eofq = eofq * -1
                                                             
    eofn.append(eofq)
    pcn.append(pcq)
    
eoff = np.asarray(eofn)
pcc = np.asarray(pcn)

############################################################################
############################################################################
############################################################################
#
#eof_wn = []
#eof_spn = []
#eof_sun = []
#eof_fn = []
#eofclimo_wn = []
#eofclimo_spn = []
#eofclimo_sun = []
#eofclimo_fn = []
#pc_wn = []
#pc_spn = []
#pc_sun = []
#pc_fn = []
#pcclimo_wn = []
#pcclimo_spn = []
#pcclimo_sun = []
#pcclimo_fn = []
#for i in xrange(slp.shape[0]):
#    eof_wq,pc_wq = calcSeasonalEOF(anomslp[i],yearslens,year1,year2,
#                    np.asarray([0,1,2]),2,2)
#    eofclimo_wq,pcclimo_wq = calcSeasonalEOF(anomslp[i],yearslens,1948,2015,
#                    np.asarray([0,1,2]),2,2)
#    eof_spq,pc_spq = calcSeasonalEOF(anomslp[i],yearslens,year1,year2,
#                    np.asarray([3,4,5]),2,2)
#    eofclimo_spq,pcclimo_spq = calcSeasonalEOF(anomslp[i],yearslens,1948,2015,
#                    np.asarray([3,4,5]),2,2)
#    eof_suq,pc_suq = calcSeasonalEOF(anomslp[i],yearslens,year1,year2,
#                    np.asarray([6,7,8]),2,2)
#    eofclimo_suq,pcclimo_suq = calcSeasonalEOF(anomslp[i],yearslens,1948,2015,
#                    np.asarray([6,7,8]),2,2)
#    eof_fq,pc_fq = calcSeasonalEOF(anomslp[i],yearslens,year1,year2,
#                    np.asarray([9,10,11]),2,2)
#    eofclimo_fq,pcclimo_fq = calcSeasonalEOF(anomslp[i],yearslens,1948,2015,
#                    np.asarray([9,10,11]),2,2)  
#                    
#    eof_wn.append(eof_wq)
#    eof_spn.append(eof_spq)
#    eof_sun.append(eof_suq)
#    eof_fn.append(eof_fq)
#    eofclimo_wn.append(eofclimo_wq)
#    eofclimo_spn.append(eofclimo_spq)
#    eofclimo_sun.append(eofclimo_suq)
#    eofclimo_fn.append(eofclimo_fq) 
#
#    pc_wn.append(pc_wq)
#    pc_spn.append(pc_spq)
#    pc_sun.append(pc_suq)
#    pc_fn.append(pc_fq)
#    pcclimo_wn.append(pcclimo_wq)
#    pcclimo_spn.append(pcclimo_spq)
#    pcclimo_sun.append(pcclimo_suq)
#    pcclimo_fn.append(pcclimo_fq) 
#
#    eof_w = np.asarray(eof_wn)
#    eof_sp = np.asarray(eof_spn)
#    eof_su = np.asarray(eof_sun)
#    eof_f = np.asarray(eof_fn)
#    eofclimo_w = np.asarray(eofclimo_wn)
#    eofclimo_sp = np.asarray(eofclimo_spn)
#    eofclimo_su = np.asarray(eofclimo_sun)
#    eofclimo_f = np.asarray(eofclimo_fn)
#
#    pc_w = np.asarray(pc_wn)
#    pc_sp = np.asarray(pc_spn)
#    pc_su = np.asarray(pc_sun)
#    pc_f = np.asarray(pc_fn)
#    pcclimo_w = np.asarray(pcclimo_wn)
#    pcclimo_sp = np.asarray(pcclimo_spn)
#    pcclimo_su = np.asarray(pcclimo_sun)
#    pcclimo_f = np.asarray(pcclimo_fn)     
#    
############################################################################
############################################################################
############################################################################
#    
#AOindex_w = (pc_w[:,:,0] - np.nanmean(pcclimo_w[:,:,0]))/np.std(pcclimo_w[:,:,0])    
#AOindex_sp = (pc_sp[:,:,0] - np.nanmean(pcclimo_sp[:,:,0]))/np.std(pcclimo_sp[:,:,0]) 
#AOindex_su = (pc_su[:,:,0] - np.nanmean(pcclimo_su[:,:,0]))/np.std(pcclimo_su[:,:,0]) 
#AOindex_f = (pc_f[:,:,0] - np.nanmean(pcclimo_f[:,:,0]))/np.std(pcclimo_f[:,:,0])  
#
############################################################################
############################################################################
############################################################################
#
#### Calculate Dipole
#def ArcticDipoleIndex(anomslp,eofpattern,years,year1,year2,monthind):
#    """
#    Calculate Arctic Dipole index using Overland et al. [2010]
#    convention by regressing slp anomalies to EOF spatial pattern
#    and then averaging per season
#    """
#    print '\n>>> Using Arctic DipoleIndex function!'       
#
#    ### Slice years
#    if np.isfinite(year1):
#        if np.isfinite(year2):
#            yearqq = np.where((years >= year1) & (years <= year2))
#            anomslp = anomslp[yearqq,:,:,:].squeeze()
#        else:
#            print 'Using entire time series for this EOF!'
#    else:
#        print 'Using entire time series for this EOF!'   
#    print 'Sliced time period for seasonal mean!'  
#    
#    ad = np.empty((anomslp.shape[0],anomslp.shape[1]))
#    for i in xrange(anomslp.shape[0]):
#        print 'Regressing year ---> %s!' % (np.arange(year1,year2+1,1)[i])
#        for j in xrange(anomslp.shape[1]):
#            varx = np.ravel(anomslp[i,j,:,:])
#            vary = np.ravel(eofpattern[1,:,:])
#            mask = np.isfinite(varx) & np.isfinite(vary)     
#            
#            ad[i,j],intercept,r,p_value,std_err = sts.stats.linregress(varx[mask],
#                                                                  vary[mask]) 
#                                                                  
#    ### Average over months
#    ad = ad[:,monthind]
#    ad = np.nanmean(ad[:,:],axis=1)     
#                                                         
#    print '*Completed: AD Regressed!\n'                                                             
#    return ad
#    
#    
#year1a = 1920
#year2a = 2080
#    
#ad_wn = []
#ad_spn = []
#ad_sun = []
#ad_fn = []
#adclimo_wn = []  
#adclimo_spn = []
#adclimo_sun = []
#adclimo_fn = []
#for i in xrange(slp.shape[0]):
#    ad_wq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,year1a,
#                       year2a,np.asarray([0,1,2]))  
#    adclimo_wq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,1981,
#                           2010,np.asarray([0,1,2]))   
#    ad_spq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,year1a,
#                           year2a,np.asarray([3,4,5]))  
#    adclimo_spq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,1981,
#                           2010,np.asarray([3,4,5]))  
#    ad_suq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,year1a,
#                           year2a,np.asarray([6,7,8]))  
#    adclimo_suq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,1981,
#                           2010,np.asarray([6,7,8]))  
#    ad_fq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,year1a,
#                           year2a,np.asarray([9,10,11]))  
#    adclimo_fq = ArcticDipoleIndex(anomslp[i],eoff[i],yearslens,1981,
#                           2010,np.asarray([9,10,11])) 
#                           
#    ad_wn.append(ad_wq)
#    ad_spn.append(ad_spq)
#    ad_sun.append(ad_suq)
#    ad_fn.append(ad_fq)
#    adclimo_wn.append(adclimo_wq)
#    adclimo_spn.append(adclimo_spq)
#    adclimo_sun.append(adclimo_suq)
#    adclimo_fn.append(adclimo_fq) 
#
#ad_w = np.asarray(ad_wn)
#ad_sp = np.asarray(ad_spn)
#ad_su = np.asarray(ad_sun)
#ad_f = np.asarray(ad_fn) 
#adclimo_w = np.asarray(adclimo_wn)
#adclimo_sp = np.asarray(adclimo_spn)
#adclimo_su = np.asarray(adclimo_sun)
#adclimo_f = np.asarray(adclimo_fn)   
#
#adindex_w = (ad_w - np.nanmean(adclimo_w))/np.std(adclimo_w)  
#adindex_sp = (ad_sp - np.nanmean(adclimo_sp))/np.std(adclimo_sp)  
#adindex_su = (ad_su - np.nanmean(adclimo_su))/np.std(adclimo_su)  
#adindex_f = (ad_f - np.nanmean(adclimo_f))/np.std(adclimo_f)  

############################################################################
############################################################################ 
############################################################################            
#### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for i in xrange(eoff.shape[0]):
    
    varf = eoff[i,0,:,:]    
    
    ax = plt.subplot(7,6,i+1)
    
    m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
                resolution='l',round =True)
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.3)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0,color='k',fontsize=6)
    m.drawmeridians(meridians,labels=[False,False,False,False],
                    linewidth=0,color='k',fontsize=6)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    # Make the plot continuous
    barlim = np.arange(-3,4,1)
    values = np.arange(-3,3.1,0.25)
    
    varf, lons_cyclic = addcyclic(varf, lons)
    varf, lons_cyclic = shiftgrid(180., varf, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
    x, y = m(lon2d, lat2d)
    
    cs = m.contourf(x,y,varf,
                    values,extend='both')
    cs1 = m.contour(x,y,varf,
                    values,linewidths=0.1,colors='darkgrey',
                    linestyles='-')
                    
    ax.annotate(r'\textbf{%s}' % ense[i], xy=(0, 0), xytext=(-0.3, 0.9),
            xycoords='axes fraction',fontsize=8,color='darkgrey',
            rotation=0)
            
    cmap = ncm.cmap('nrl_sirkes')         
    cs.set_cmap(cmap) 

cbar_ax = fig.add_axes([0.51,0.16,0.31,0.027])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))  
cbar.ax.tick_params(labelsize=8)   
cbar.set_label(r'\textbf{hPa}')

plt.subplots_adjust(wspace=-0.6)

fig.suptitle(r'\textbf{EOF1, NDJFM 1981-2010 (SLP) -- LENS}')

plt.savefig(directoryfigure + 'testeof1.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################

#fig = plt.figure()
#ax = plt.subplot(221)
#
#### Plot PC time series
#### Adjust axes in time series plots 
#def adjust_spines(ax, spines):
#    for loc, spine in ax.spines.items():
#        if loc in spines:
#            spine.set_position(('outward', 10))
#        else:
#            spine.set_color('none')  
#    if 'left' in spines:
#        ax.yaxis.set_ticks_position('left')
#    else:
#        ax.yaxis.set_ticks([])
#
#    if 'bottom' in spines:
#        ax.xaxis.set_ticks_position('bottom')
#    else:
#        ax.xaxis.set_ticks([])
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#zeros = [0] * adindex_w.shape[1]
#time = np.arange(adindex_w.shape[1])
#
#plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=3)
#
#plt.plot(time,adindex_w.transpose(),zorder=2,linewidth=0.5,
#         color='steelblue',label='PC2',alpha=0.2)
#         
#plt.plot(time,np.nanmean(adindex_w,axis=0),zorder=4,linewidth=1.3,
#         color='darkred',label='PC2',alpha=1)
#
#plt.axvline(85,linestyle='--',linewidth=1.2,color='k') 
#         
#plt.xticks(np.arange(0,year2a-year1a+3,20),
#           map(str,np.arange(year1a,year2a+3,20)),fontsize=8)
#plt.xlim([0,year2a-year1a])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{JFM}',fontsize=20,color='darkgrey')
#
#ax = plt.subplot(222)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#zeros = [0] * adindex_sp.shape[1]
#
#plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=3)
#
#plt.plot(time,adindex_sp.transpose(),zorder=2,linewidth=0.5,
#         color='steelblue',label='PC2',alpha=0.2)
#         
#plt.plot(time,np.nanmean(adindex_sp,axis=0),zorder=4,linewidth=1.3,
#         color='darkred',label='PC2',alpha=1)
#         
#plt.axvline(85,linestyle='--',linewidth=1.2,color='k')         
#         
#plt.xticks(np.arange(0,year2a-year1a+3,20),
#           map(str,np.arange(year1a,year2a+3,20)),fontsize=8)
#plt.xlim([0,year2a-year1a])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{AMJ}',fontsize=20,color='darkgrey')
#
#
#ax = plt.subplot(223)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#zeros = [0] * adindex_su.shape[1]
#
#plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=3)
#
#plt.plot(time,adindex_su.transpose(),zorder=2,linewidth=0.5,
#         color='steelblue',label='PC2',alpha=0.2)
#         
#plt.plot(time,np.nanmean(adindex_su,axis=0),zorder=4,linewidth=1.3,
#         color='darkred',label='PC2',alpha=1)
#
#plt.axvline(85,linestyle='--',linewidth=1.2,color='k') 
#         
#plt.xticks(np.arange(0,year2a-year1a+3,20),
#           map(str,np.arange(year1a,year2a+3,20)),fontsize=8)
#plt.xlim([0,year2a-year1a])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{JAS}',fontsize=20,color='darkgrey')
#
#
#ax = plt.subplot(224)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
#zeros = [0] * adindex_f.shape[1]
#
#plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=3)
#
#plt.plot(time,adindex_f.transpose(),zorder=2,linewidth=0.5,
#         color='steelblue',label='PC2',alpha=0.2)
#         
#plt.plot(time,np.nanmean(adindex_f,axis=0),zorder=4,linewidth=1.3,
#         color='darkred',label='PC2',alpha=1)
#
#plt.axvline(85,linestyle='--',linewidth=1.2,color='k') 
#         
#plt.xticks(np.arange(0,year2a-year1a+3,20),
#           map(str,np.arange(year1a,year2a+3,20)),fontsize=8)
#plt.xlim([0,year2a-year1a])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{OND}',fontsize=20,c############################################################################
############################################################################ 
############################################################################            
#### Plot figure
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()

for i in xrange(eoff.shape[0]):
    
    varf = eoff[i,0,:,:]    
    
    ax = plt.subplot(7,6,i+1)
    
    m = Basemap(projection='npstere',boundinglat=70,lon_0=270,
                resolution='l',round =True)
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='k',linewidth=0.3)
    parallels = np.arange(50,90,10)
    meridians = np.arange(-180,180,30)
    m.drawparallels(parallels,labels=[False,False,False,False],
                    linewidth=0,color='k',fontsize=6)
    m.drawmeridians(meridians,labels=[False,False,False,False],
                    linewidth=0,color='k',fontsize=6)
    m.drawlsmask(land_color='darkgrey',ocean_color='mintcream')
    
    # Make the plot continuous
    barlim = np.arange(-3,4,1)
    values = np.arange(-3,3.1,0.25)
    
    varf, lons_cyclic = addcyclic(varf, lons)
    varf, lons_cyclic = shiftgrid(180., varf, lons_cyclic, start=False)
    lon2d, lat2d = np.meshgrid(lons_cyclic, lats)
    x, y = m(lon2d, lat2d)
    
    cs = m.contourf(x,y,varf,
                    values,extend='both')
    cs1 = m.contour(x,y,varf,
                    values,linewidths=0.1,colors='darkgrey',
                    linestyles='-')
                    
    ax.annotate(r'\textbf{%s}' % ense[i], xy=(0, 0), xytext=(-0.3, 0.9),
            xycoords='axes fraction',fontsize=8,color='darkgrey',
            rotation=0)
            
    cmap = ncm.cmap('nrl_sirkes')         
    cs.set_cmap(cmap) 

cbar_ax = fig.add_axes([0.51,0.16,0.31,0.027])                
cbar = fig.colorbar(cs,cax=cbar_ax,orientation='horizontal',
                    extend='both',extendfrac=0.07,drawedges=True)                      
cbar.set_ticks(barlim)
cbar.set_ticklabels(map(str,barlim))  
cbar.ax.tick_params(labelsize=8)   
cbar.set_label(r'\textbf{hPa}')

plt.subplots_adjust(wspace=-0.6)

fig.suptitle(r'\textbf{EOF1, NDJFM 1981-2010 (SLP) -- LENS}')

plt.savefig(directoryfigure + 'testeof1.png',dpi=300)
#
#
#fig.subplots_adjust(hspace=0.4)
#plt.savefig(directoryfigure + 'adindex_LENS_all.png',dpi=300)
#
############################################################################
############################################################################
############################################################################
#### Create text files
#directorytext = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'
#
#np.savetxt(directorytext + 'AD_JFM_%s%s_LENS.txt' % (year1a,year2a),
#           adindex_w)
#np.savetxt(directorytext + 'AD_AMJ_%s%s_LENS.txt' % (year1a,year2a),
#           adindex_sp)
#np.savetxt(directorytext + 'AD_JAS_%s%s_LENS.txt' % (year1a,year2a),
#           adindex_su)
#np.savetxt(directorytext + 'AD_OND_%s%s_LENS.txt' % (year1a,year2a),
#           adindex_f)  
#         
#plt.bar(np.arange(len(AOindex_su))-0.4,AOindex_su,label='PC1',color='indianred')
#
#plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
#plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
#         
#plt.xticks(np.arange(0,year2-year1+3,5),
#           map(str,np.arange(year1,year2+3,5)),fontsize=6)
#plt.xlim([0,year2-year1+2])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{JAS}',fontsize=20)
#
#
#ax = plt.subplot(224)
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#plt.grid(color='w',zorder=1,alpha=0.2)
#
#zeros = [0] * len(AOindex_f)
#
#plt.plot(zeros,linewidth=1.5,color='k',linestyle='-',zorder=1)
#
#plt.plot(adindex_f,marker='o',markersize=4,zorder=2,linewidth=2,
#         color='steelblue',label='PC2')
#         
#plt.bar(np.arange(len(AOindex_f))-0.4,AOindex_f,label='PC1',color='indianred')
#
#plt.axvline(np.where(years2 == 2007)[0],color='k',linestyle='--')
#plt.axvline(np.where(years2 == 2012)[0],color='k',linestyle='--')
#         
#plt.xticks(np.arange(0,year2-year1+3,5),
#           map(str,np.arange(year1,year2+3,5)),fontsize=6)
#plt.xlim([0,year2-year1+2])
#
#plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=6)
#plt.ylim([-3,3])
#
#plt.text(-1.3,2.85,r'\textbf{OND}',fontsize=20)
#
#
#fig.subplots_adjust(hspace=0.4)
#plt.savefig(directoryfigure + 'eof2pattern_8110.png',dpi=300)