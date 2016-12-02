"""
Script looks at NCEP/NCAR reanalysis trends of SLP with SIV/SIT
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 27 October 2016
"""

### Import modules
import numpy as np
import matplotlib.pyplot as plt
import datetime
import scipy.stats as sts
from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid
import read_SeaIceThick_PIOMAS as CT
import read_SeaIceConc_PIOMAS as CC
import calc_PiomasArea as CA
import calc_SIV as CV
import read_NCEP as NP

### Define directories
directorydata = '/home/zlabe/Surtsey/NCEP/'
directorydata2 =   
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
yearmin = 1979
yearmax = 2015
years = np.arange(yearmin,yearmax+1,1)
months = [r'Jan',r'Feb',r'Mar',r'Apr',r'May',r'Jun',r'Jul',r'Aug',
          r'Sep',r'Oct',r'Nov',r'Dec']
          
### Read in functions
lats,lons,slp = NP.readNCEP(directorydata,years,'slp','surface')
lats,lons,sit = CT.readPiomas(directorydata,years,0.15)
lats,lons,sic = CC.readPiomas(directorydata,years,0.01)
area = CA.readPiomasArea(directorydata)
sivq = CV.sivGrid(sit,sic,area,False)    