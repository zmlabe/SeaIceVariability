"""
Function reads NCEP/NCAR reanalysis for a variety of variables according
to the reader for monthly data.
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 28 September 2016
    
Usage
-----
    readNCEPWind(directory,years,wind)
"""

def readNCEPWind(directory,years,wind):
    """
    Function reads monthly NCEP/NCAR u/v wind reanalysis data.

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    years : integers
        years for data files
    wind: string
        u,v

    Returns
    -------
    lats : 2d array
        latitudes
    lons : 2d array
        longitudes
    var : 4d array [year,month,lat,lon]
        ncep variable

    Usage
    -----
    lats,lons,var = readNCEP(directory,years,name,wind)
    """
    
    print '\n>>> Using readNCEP Wind function!'    
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
    if wind == 'u':
        directoryname = directory + 'wind/'
        filename = directoryname + 'uwnd.mon.mean.nc'
        varname = 'uwnd'
    elif wind == 'v':
        directoryname = directory + 'wind/'
        filename = directoryname + 'vwnd.mon.mean.nc'
        varname = 'vwnd'
    else:
        ValueError('Wrong name of variable!')
        
    print 'Finished reading %s wind data!' % wind
        
    ### Time dimension
    yearstart = 1948                               # January 1
    yearmin = int(((np.nanmin(years)-yearstart)*12.))
    yearmax = int((np.nanmax(years)-yearstart)*12.)+12
    print 'Slicing years %s-%s!' % (np.nanmin(years),np.nanmax(years))
        
    ### Read netcdf file
    data = Dataset(filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    varq = data.variables[varname][yearmin:yearmax,:,:]
    data.close()   
    print 'Finished reading data file!'
    
    ### Reshaping array
    var = np.reshape(varq,(varq.shape[0]/12.,12,
                           lats.shape[0],lons.shape[0]))
    
    
    print '*Completed: NCEP wind data read!\n'
    return lats,lons,var
