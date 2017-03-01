"""
Function reads 20CRv2 reanalysis for a variety of variables according
to the reader for monthly data. Data is available from 01/1851
through 12/2014
 
Notes
-----
    Source : http://www.esrl.noaa.gov/psd/data/gridded/data.
             ncep.reanalysis.derived.html
    Author : Zachary Labe
    Date   : 17 February 2017
    
Usage
-----
    read20CR(directory,years,name,height)
"""

def read20CR(directory,years,name,height):
    """
    Function reads monthly NCEP/NCAR reanalysis data.

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    years : integers
        years for data files
    name : string
        prmsl
        
    Returns
    -------
    lats : 2d array
        latitudes
    lons : 2d array
        longitudes
    var : 4d array [year,month,lat,lon]
        20CR variable

    Usage
    -----
    lats,lons,var = read20CR(directory,years,name,height)
    """
    
    print '\n>>> Using readNCEP function!'    
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
            
    if name == 'prmsl':
        directoryname = directory + 'prmsl/'
        filename = directoryname + 'prmsl.mon.mean.nc'
        varname = 'prmsl'
        if height != 'surface':
            ValueError('Wrong height value given!')
    else:
        ValueError('Wrong name of variable!')
        
    ### Time dimension
    yearstart = 1851                               # January 1
    yearmin = int(((np.nanmin(years)-yearstart)*12.))
    yearmax = int((np.nanmax(years)-yearstart)*12.)+12
    print 'Slicing years %s-%s!' % (np.nanmin(years),np.nanmax(years))
        
    ### Read netcdf file
    data = Dataset(filename)
    lats = data.variables['lat'][:]
    lons = data.variables['lon'][:]
    if name == 'prmsl':
        varq = data.variables[varname][yearmin:yearmax,:,:]
        
        ### Conversion Pa to hPa
        varq = varq/100.
    else:
        varq = data.variables[varname][yearmin:yearmax,z,:,:]
    data.close()   
    print 'Finished reading data file!'
    
    ### Reshaping array
    var = np.reshape(varq,(varq.shape[0]/12.,12,
                           lats.shape[0],lons.shape[0]))
        
    print '*Completed: NCEP data read!\n'
    return lats,lons,var