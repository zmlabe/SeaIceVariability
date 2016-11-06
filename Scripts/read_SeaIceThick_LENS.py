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
    Function reads LENS control netCDF4 data array

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
    lats,lons,sit = readLENS(directory,years,threshold)
    """
    
    print '\n>>> Using readLENS function!'
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ### Modify directory
    directory = directory + 'CESM_large_ensemble/SIT/interp_1deg/'
    files = 'SIT_control_0001-1800.nc'
    filename = directory + files
    
    data = Dataset(filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    sitq = data.variables['SIT'][:,:,:]
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
    
def readLENSEnsemble(directory,threshold):
    """
    Function reads LENS ensembles netCDF4 data array

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
    sit : 5d array [ens,year,month,lat,lon]
        sea ice thickness (m) 

    Usage
    -----
    lats,lons,sit = readLENS(directory,years,threshold)
    """
    
    print '\n>>> Using readLENS function!'
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    ens = np.array(['02','03','04','05','06','07','08','09'])
#    ens = np.array(['02','03','04','05','06','07','08','09']) + \
#        map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))
    
    ### Modify directory
    directory = directory + 'CESM_large_ensemble/SIT/interp_1deg/'
    
    sitq = np.empty((len(ens),86*12,24,360))
    for i in xrange(len(ens)):
        files = 'b.e11.B20TRC5CNBDRD.f09_g16.0%s.cice.h.hi_nh.192001-200512.nc' % ens[i]
        filename = directory + files
        
        data = Dataset(filename)
        lats = data.variables['lat'][156:180]
        lons = data.variables['lon'][:]
        sitq[i,:,:,:] = data.variables['SIT'][:,156:180,:] # lats > 65
        data.close()
        
        print 'Completed reading LENS ensemble #%s!' % ens[i]
        
    sit = np.reshape(sitq,(len(ens),sitq.shape[1]/12,12,
                           lats.shape[0],lons.shape[0]))
    sit = np.squeeze(np.asarray(sit))
    
    ### Mask out threshold values
    if threshold == 'None':
        sit[np.where(sit < 0)] = np.nan
        sit[np.where(sit > 12)] = np.nan
    else:
        sit[np.where(sit < threshold)] = np.nan
        sit[np.where(sit < 0)] = np.nan
        sit[np.where(sit > 12)] = np.nan
        
    print 'Masking SIT data < %s m!' % threshold

    print '*Completed: Read SIT data!' 
    
    return sit,lats,lons
    
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