# SeaIceVariability [![DOI](https://zenodo.org/badge/68020882.svg)](https://zenodo.org/badge/latestdoi/68020882)
Arctic sea ice interannual variability and change

###### Under construction... ```[Python 2.7]```

## Contact
Zachary Labe - [Research Website](http://sites.uci.edu/zlabe/) - [@ZLabe](https://twitter.com/ZLabe)


## Description
Changes in Arctic sea ice extent and concentration are well documented within the satellite record (1979). However, quantifying pan-Arctic sea ice thickness and volume is challenging as a result of limited observations until recently through ICESat (2003-2009) and CryoSat-2 (2011-present). The Pan-Arctic Ice-Ocean Modeling and Assimilation System (PIOMAS) is a POP ocean and sea ice model that has been widely validated in reproducing a sea ice thickness and volume record consistent with the spatial and temporal variability of our limited observations. This project will provide a comprehensive overview of long-term trends and variability in sea ice volume using PIOMAS in addition to ICESat/CryoSat-2. We are particularly interested in the the distribution of sea ice and regional trends, which may be important for feedbacks with the large-scale atmospheric circulation.

+ ```Data/```: Additional data files (ascii and netCDF4) modified from original sources. Data includes calculations from EOF analysis, linear trends, and nearest-neighbor interpolation 
+ ```Scripts/```: Main [Python](https://www.python.org/) scripts/functions used in data analysis and plotting. These scripts are not compatible with Python 3+.
+ ```requirements.txt```: List of environments and modules associated with the most recent version of this project. A Python [Anaconda2 Distribution](https://docs.continuum.io/anaconda/) was used for our analysis. All AGCM experiments were processed through resources on CISL's [Cheyenne](https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne) supercomputer. Tools including [NCL](https://www.ncl.ucar.edu/), [CDO](https://code.mpimet.mpg.de/projects/cdo), and [NCO](http://nco.sourceforge.net/) were also used for initial data manipulation.

## Data
+ CESM Large Ensemble Project (LENS) : [[DATA]](http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html)
    + Kay, J. E and Coauthors, 2015: The Community Earth System Model (CESM) large ensemble project: A community resource for studying climate change in the presence of internal climate variability. Bull. Amer. Meteor. Soc., 96, 1333–1349, doi:10.1175/BAMS-D-13-00255.1, [[Publication]](http://journals.ametsoc.org/doi/full/10.1175/BAMS-D-13-00255.1)
+ CryoSat-2 : [[DATA]](http://data.meereisportal.de/gallery/index_new.php?ice-type=thickness&satellite=C&region=n&resolution=monthly&minYear=2010&minMonth=11&maxYear=2010&maxMonth=11&showMaps=y&dateRepeat=n&submit2=Show+&lang=de_DE&active-tab=thickness)
    + Laxon, S. W., and Coauthors, 2013: CryoSat-2 estimates of Arctic sea ice thickness and volume. Geophysical Research Letters, 40 (4), 732–737, doi:10.1002/grl.50193, [[Publication]](http://doi.wiley.com/10.1002/grl.50193)
+ Ice, Cloud, and land Elevation Satellite (ICESat) : [[DATA]](https://rkwok.jpl.nasa.gov/icesat/)
    + Kwok, R., 2004: ICESat observations of Arctic sea ice: A first look. Geophysical Research Letters, 31 (16), L16 401, doi:10.1029/2004GL020309, [[Publication]](http://doi.wiley.com/10.1029/2004GL020309)
+ Pan-Arctic Ice Ocean Modeling and Assimilation System (PIOMAS) : [[DATA]](http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid)
    + Zhang, J., and D. A. Rothrock, 2003: Modeling Global Sea Ice with a Thickness and Enthalpy Distribution Model in Generalized Curvilinear Coordinates. Monthly Weather Review, 131 (5), 845–861, doi:10.1175/1520-0493(2003)131<0845:MGSIWA>2.0.CO;2 [[Publication]](http://journals.ametsoc.org/doi/abs/10.1175/1520-0493%282003%29131%3C0845%3AMGSIWA%3E2.0.CO%3B2)
+ Submarine sea ice thickness : [[DATA]](http://psc.apl.uw.edu/sea_ice_cdr/)
    + Lindsay, R., 2010: New Unified Sea Ice Thickness Climate Data Record. Eos, Transactions Amer578
ican Geophysical Union, 91 (44), 405, doi:10.1029/2010EO440001, [[Publication]](http://doi.wiley.com/10.1029/2010EO440001)

## Publications
+ **Labe, Z.M.**, G. Magnusdottir, and H.S. Stern, 2018. Variability of Arctic sea ice thickness using PIOMAS and the CESM Large Ensemble, *in revision*

## Conferences
+ [2] **Labe, Z.M.**, G. Magnusdottir, and H.S. Stern. Variability and future projections of Arctic sea ice thickness. *Understanding the Causes and Consequences of Polar Amplification Workshop*, Aspen Global Change Institute, Aspen, CO (Jun 2017). [[RECORDING]](https://www.agci.org/lib/17s1/variability-and-future-projections-arctic-sea-ice-thickness) [[SLIDES]](https://www.agci.org/lib/17s1/variability-and-future-projections-arctic-sea-ice-thickness)
+ [1] **Labe, Z.M.**, G. Magnusdottir, and H.S. Stern. SArctic Sea Ice Thickness Variability and the Large-scale Atmospheric Circulation Using Satellite Observations, PIOMAS, and the CESM Large Ensemble, 14th Conference on Polar Meteorology and Oceanography, 97th Annual Meeting of the American Meteorological Society (Jan 2017). [[ABSTRACT]](https://ams.confex.com/ams/97Annual/webprogram/Paper313445.html) [[POSTER]](http://sites.uci.edu/zlabe/files/2017/01/ZLabeMagnusdottirStern_AMSposter_2017.pdf)
