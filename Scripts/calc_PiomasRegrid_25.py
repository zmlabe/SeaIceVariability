"""
Function regrids Piomas SIT data onto 25km EASE grid.
 
Notes
-----
    Source : http://psc.apl.washington.edu/zhang/IDAO/data_piomas.html
    Author : Zachary Labe
    Date   : 13 September 2016
    
Usage
-----
    sitnew,lats,lons = regridPiomas(directory,sit,lat1,lon1,lat2,lon2)
"""

print '\n>>> Using regridPiomas function!'

import numpy as np
import calc_grid25 as C25
import read_SeaIceThick_PIOMAS as CT

### Define directories
directorydata1 = '/home/zlabe/Surtsey/seaice_obs/PIOMAS/' 
directorydata2 = '/home/zlabe/Surtsey/seaice_obs/sic/sic_rename/' 

### Alott time series
yearmin = 1979
yearmax = 2016
years = np.arange(yearmin,yearmax+1,1)

### Call functions
lat2,lon2 = C25.readGrid25(directorydata2)
lat1,lon1,sit = CT.readPiomas(directorydata1,years,0.0)