"""
Scripts plots Arctic Dipole for LENS (1920-2080) averaged by 3-month
seasons
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Reference : Wu et al. [2006] and Overland et al. [2012]
    Author : Zachary Labe
    Date   : 21 February 2017
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap
import statsmodels.api as sm

### Define directories
directorydata = '/home/zlabe/Documents/Research/SeaIceVariability/Data/'  
directoryfigure = '/home/zlabe/Desktop/'

### Define time           
now = datetime.datetime.now()
currentmn = str(now.month)
currentdy = str(now.day)
currentyr = str(now.year)
currenttime = currentmn + '_' + currentdy + '_' + currentyr
titletime = currentmn + '/' + currentdy + '/' + currentyr
print '\n' '----Plot LENS AD - %s----' % titletime 

### Alott time series
year1 = 1920
year2 = 2080
years = np.arange(year1,year2+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in data
seasons = ['JFM','AMJ','JAS','OND']  

ad = np.empty((4,39,len(years)))
for i in xrange(len(seasons)):
    filename1 = 'AD_%s_%s%s_LENS.txt' % (seasons[i],year1,year2)
    ad[i,:,:] = np.genfromtxt(directorydata + filename1)

    print 'Reading data file %s!' % seasons[i]

timex = np.arange(len(ad[0,0,:]))
vary1 = ad[0,0,:] 
slope1,intercept1,r1,p_value,std_err = sts.stats.linregress(timex,vary1) 
line1 = slope1*timex + intercept1
#vary2 = ad[1,:] 
#slope2,intercept2,r2,p_value,std_err = sts.stats.linregress(timex,vary2) 
#line2 = slope2*timex + intercept2
#vary3 = ad[2,:] 
#slope3,intercept3,r3,p_value,std_err = sts.stats.linregress(timex,vary3)
#line3 = slope3*timex + intercept3 
#vary4 = ad[3,:] 
#slope4,intercept4,r4,p_value,std_err = sts.stats.linregress(timex,vary4) 
#line4 = slope4*timex + intercept4

smoothed1 = sm.nonparametric.lowess(vary1,timex)
#smoothed2 = sm.nonparametric.lowess(vary2,timex)
#smoothed3 = sm.nonparametric.lowess(vary3,timex)
#smoothed4 = sm.nonparametric.lowess(vary4,timex)

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
        
plt.rc('text',usetex=True)
plt.rc('font',**{'family':'sans-serif','sans-serif':['Avant Garde']}) 

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

zeros = [0] * ad.shape[2]

plt.plot(zeros,linewidth=1,color='k',linestyle='-',zorder=1)

plt.bar(np.arange(ad.shape[2])-0.4,ad[0,0,:],color='lightskyblue',edgecolor='lightskyblue')
plt.plot(timex,line1,zorder=1,linewidth=2,color='indianred')
plt.plot(smoothed1[:,0],smoothed1[:,1],color='darkslateblue',zorder=8,linewidth=1.2)
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,year2-year1+3,10),
           map(str,np.arange(year1,year2+3,10)),fontsize=8)
plt.xlim([0,year2-year1])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
plt.ylim([-3,3])

plt.text(-1.3,2.85,r'\textbf{JAS}',fontsize=20,color='darkgrey')


plt.savefig(directoryfigure + 'ad_LENS_test.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('darkgrey')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')

zeros = [0] * ad.shape[2]

plt.plot(ad[3,0,:].transpose(),color='darkorchid',linewidth=0.3)
plt.plot(ad[3,1,:].transpose(),color='darkorange',linewidth=0.3)
plt.plot(ad[3,2,:].transpose(),color='indianred',linewidth=0.3)
plt.plot(ad[3,3,:].transpose(),color='darkgreen',linewidth=0.3)
plt.plot(ad[3,4,:].transpose(),color='cornflowerblue',linewidth=0.3)

plt.plot(np.nanmean(ad[2,:5,:],axis=0),color='r',linewidth=2)


plt.plot(zeros,linewidth=1,color='k',linestyle='-',zorder=2)
plt.axvline(85,linestyle='--',linewidth=2,color='k')

plt.xticks(np.arange(0,year2-year1+3,10),
           map(str,np.arange(year1,year2+3,10)),fontsize=8)
plt.xlim([0,year2-year1])

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)),fontsize=8)
plt.ylim([-3,3])

plt.savefig(directoryfigure + 'ad_LENS_seasontest.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
#fig = plt.figure()
#ax = plt.subplot(111)
#
#test = ad[3,1,:]
#
#### Adjust axes spines
#adjust_spines(ax, ['left', 'bottom'])
#ax.spines['top'].set_color('none')
#ax.spines['right'].set_color('none')
#ax.spines['left'].set_color('darkgrey')
#ax.spines['bottom'].set_color('darkgrey')
#ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
#
##n, bins, patches = ax.hist(sts.norm(ad[3,:,:],bins=range(-3,3),align='left')
#fit = sts.norm.pdf(test,np.nanmean(test),np.nanstd(test))
#
#plt.plot(test,fit,'-y')
#plt.hist(test,normed=True)
#
##for i in range(len(patches)):
##    patches[i].set_facecolor('mediumseagreen')
##    patches[i].set_edgecolor('white')
##    patches[i].set_linewidth(0.9)
#
##plt.xticks(np.arange(-1,1.1,0.25),np.arange(-1,1.1,0.25))
##plt.yticks(np.arange(0,7,1),map(str,np.arange(0,7,1)))
##plt.xlim([-1,1])
##plt.ylim([0,6])
#
##plt.ylabel(r'\textbf{Number of Occurences}')
#
#plt.savefig(directoryfigure + 'test_ad_hist.png',dpi=300)

###########################################################################
###########################################################################
###########################################################################
adq = np.reshape(ad[:,:,11:],(ad.shape[0],ad.shape[1]*25.,ad.shape[2]/25.))\
#adq = np.reshape(ad[:,:,11:],(ad.shape[0],ad.shape[1],25,ad.shape[2]/25.))

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('white')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params(axis='x',bottom='off',labelbottom='off') 

bx = plt.boxplot(adq[3,:,:],0,'',patch_artist=True)

plt.text(0.7,-3.2,r'\textbf{1931-1955}',fontsize=8,color='darkgrey')
plt.text(1.7,-3.2,r'\textbf{1956-1980}',fontsize=8,color='darkgrey')
plt.text(2.7,-3.2,r'\textbf{1981-2005}',fontsize=8,color='darkgrey')
plt.text(3.7,-3.2,r'\textbf{2006-2030}',fontsize=8,color='darkgrey')
plt.text(4.7,-3.2,r'\textbf{2031-2055}',fontsize=8,color='darkgrey')
plt.text(5.7,-3.2,r'\textbf{2056-2080}',fontsize=8,color='darkgrey')


for i in bx['caps']:
    i.set(color='k',linewidth=0)
for whisker in bx['whiskers']:
    whisker.set(color='k',linestyle='-')
for box in bx['boxes'][:]:
    box.set(color='indianred')


plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)))    
plt.ylabel(r'\textbf{AD Index}',fontsize=13)
plt.text(0.4,2.85,r'\textbf{OND}',fontsize=20,color='darkgrey')

plt.savefig(directoryfigure + 'ad_LENS_OND_25yr.png',dpi=600)

###########################################################################
###########################################################################
###########################################################################
adq2 = np.reshape(ad[:,:,1:],(ad.shape[0],ad.shape[1]*10.,ad.shape[2]/10.))\

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('white')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params(axis='x',bottom='off',labelbottom='off') 

bx = plt.boxplot(adq2[3,:,:],0,'',patch_artist=True)

for i in bx['caps']:
    i.set(color='k',linewidth=0)
for whisker in bx['whiskers']:
    whisker.set(color='k',linestyle='-')
for box in bx['boxes'][:]:
    box.set(color='indianred')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)))    
plt.ylabel(r'\textbf{AD Index}',fontsize=13)
plt.text(0.3,2.85,r'\textbf{OND}',fontsize=20,color='darkgrey')

plt.savefig(directoryfigure + 'ad_lens_OND_10yr.png',dpi=600)

###########################################################################
###########################################################################
###########################################################################
adq3 = np.reshape(ad[:,:,1:],(ad.shape[0],ad.shape[1]*5.,ad.shape[2]/5.))\

fig = plt.figure()
ax = plt.subplot(111)

### Adjust axes spines
adjust_spines(ax, ['left', 'bottom'])
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['left'].set_color('darkgrey')
ax.spines['bottom'].set_color('white')
ax.tick_params('both',length=4,width=1.5,which='major',color='darkgrey')
ax.tick_params(axis='x',bottom='off',labelbottom='off') 

bx = plt.boxplot(adq3[3,:,:],0,'',patch_artist=True)

for i in bx['caps']:
    i.set(color='k',linewidth=0)
for whisker in bx['whiskers']:
    whisker.set(color='k',linestyle='-')
for box in bx['boxes'][:]:
    box.set(color='indianred')

plt.yticks(np.arange(-3,4,1),map(str,np.arange(-3,4,1)))    
plt.ylabel(r'\textbf{AD Index}',fontsize=13)

plt.text(0,2.85,r'\textbf{OND}',fontsize=20,color='darkgrey')

plt.savefig(directoryfigure + 'ad_lens_OND_5yr.png',dpi=600)



