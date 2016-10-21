"""
Script reads LENS data for SIT for the control run.
 
Notes
-----
    Author : Zachary Labe
    Date   : 20 October 2016
    
Usage
-----
    readLENS(directory,threshold)
    meanThick(yearmin,yearmax,years,sit)
"""

def readLENS(directory,threshold):
    """
    Function reads LENS netCDF4 data array

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    threshold : float
        mask sea ice thickness amounts < to this value

    Returns
    -------
    lats : 1d array
        latitudes
    lons : 1d array
        longitudes
    sit : 4d array [year,month,lat,lon]
        sea ice thickness (m) 

    Usage
    -----
    lats,lons,sit = readPiomas(directory,years,threshold)
    """
    
    print '\n>>> Using readPiomas function!'
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Modify directory
    directory = directory + 'CESM_large_ensemble/SIT/'
    files = 'SIT_control_0001-1800.nc'
#    files = 'b.e11.B1850C5CN.f09_g16.005.cice.h.hi_nh.040001-049912.nc'
    filename = directory + files
    
    data = Dataset(filename)
    lats = data.variables['latitude'][80:]
    lons = data.variables['longitude'][:]
    sitq = data.variables['SIT'][:,80:,:]
    data.close()
    
    ### Test for slice lats > 60
#    latq = np.where(lats > 60)[0]
    
    ### Reshape into 4d
    sit = np.reshape(sitq,(sitq.shape[0]/12,12,lats.shape[0],lons.shape[0]))
    sit = np.asarray(sit)
    
    ### Mask out threshold values
    if threshold == 'None':
        sit[np.where(sit < 0)] = np.nan
    else:
        sit[np.where(sit < threshold)] = np.nan
        sit[np.where(sit < 0)] = np.nan
    print 'Masking SIT data < %s m!' % threshold

    print '*Completed: Read SIT data!'   
    
    return lats,lons,sit
    
def meanThick(sit):
    """
    Function calculates climatological gridded average sea ice thickness

    Parameters
    ----------
    sit : 4d array [year,month,lat,lon]
        sea ice thickness

    Returns
    -------
    meansit : 3d array [month,lat,lon]
        average sit over set climatological bounds

    Usage
    -----
    meansit = meanThick(sit)
    """
    
    print '\n>>> Using climatology sit function!'
    
    ### Import modules
    import numpy as np
    
    ### Calculate average
    meansit = np.nanmean(sit,axis=0)
    
    print '*Completed: Climatology sit calculated!\n'
    return meansit