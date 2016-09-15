"""
Script reads Sea Ice Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS
Passive Microwave Data, Version 1 binary files for select variables and 
regrids according to selected grid style (e.g., NSIDC EASE grid data). 
 
Notes
-----
    Source : https://nsidc.org/data/nsidc-0051#
    Author : Zachary Labe
    Date   : 13 September 2016
    
Usage
-----
    lats,lons = readGrid25(directory)
"""

def readGrid25(directory):
    """
    Function reads sic binary data to get lat/lon

    Parameters
    ----------
    directory : string
        working directory for stored PIOMAS files

    Returns
    -------
    lats : 2d array
        latitudes
    lons : 2d array
        longitudes

    Usage
    -----
    lats,lons = readPiomas(directory,years,threshold)
    """
    
    print '\n>>> Using readGrid25 function!'
    
    ### Import modules
    import numpy as np
    
    ### Read binary lat x lon arrays
    lons = np.fromfile(directory + 'psn25lons_v3.dat',dtype='<i4')
    lons = (np.reshape(lons,(448,304)))/100000.  # Scale Factor
    lats = np.fromfile(directory + 'psn25lats_v3.dat',dtype='<i4')
    lats = (np.reshape(lats,(448,304)))/100000.  # Scale Factor
    
    print '*Completed: Read grid data!' 
    return lats,lons
