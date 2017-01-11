"""
Script reads LENS data for selected variables
 
Notes
-----
    Author : Zachary Labe
    Date   : 28 November 2016
    
Usage
-----
    lats,lons,var = readLENS(directory,varq)
"""
    
def readLENSEnsemble(directory,varq):
    """
    Function reads LENS ensembles netCDF4 data array

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files
    varq : string
        variable from LENS

    Returns
    -------
    lats : 1d array
        latitudes
    lons : 1d array
        longitudes
    varq : 5d array [ens,year,month,lat,lon]
        selected variable

    Usage
    -----
    lats,lons,var = readLENS(directory,varq)
    """
    
    print '\n>>> Using readLENS function!'
    
    ### Import modules
    import numpy as np
    from netCDF4 import Dataset
    
#    ens = ['02','03']
    ens = ['02','03','04','05','06','07','08','09'] + \
        map(str,np.arange(10,36,1)) + map(str,np.arange(101,106,1))
    
    ### Modify directory
    directory = directory + '%s/' % (varq)
    
    varn = np.empty((len(ens),161*12,11,144)) # 96 for all
    for i in xrange(len(ens)):
        filename = '%s_0%s_1920-2100.nc' % (varq,ens[i])
        
        if int(ens[i]) > 100:
            filename = '%s_%s_1920-2100.nc' % (varq,ens[i])
            
        data = Dataset(directory + filename)
        lats = data.variables['latitude'][85:] # > 70N 85 11
        lons = data.variables['longitude'][:]
        varn[i,:,:,:] = data.variables['%s' % varq][:-240,85:,:] # -2080
        data.close()
        
        print 'Completed: Read LENS Ensemble #%s - %s!' % (ens[i],varq)
                
    var = np.reshape(varn,(len(ens),varn.shape[1]/12,12,
                           lats.shape[0],lons.shape[0]))
    var = np.squeeze(np.asarray(var))
    
    ### Modify Units
    if varq == 'SLP':
        var = var/100. #Pa to hPa
        

    print '*Completed: Read %s data!' % varq
    
    return var,lats,lons
    
#var,lats,lons = readLENSEnsemble('/home/zlabe/Surtsey3/CESM_large_ensemble/' ,'SST')