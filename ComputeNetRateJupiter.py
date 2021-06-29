# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 23:45:33 2021

@author: Steven Hill
"""

def ComputeNetRateJupiter(scidata,header,IDs,positions,radii):
    from photutils import CircularAperture
    from photutils import aperture_photometry
    from photutils import CircularAnnulus
    import Meta_and_Control_Data_Operations as Meta
    from astropy.table import Table, hstack
    import pylab as pl

    #Create aperture objects    
    apertures = CircularAperture(positions, r=radii[0])
    annulus_apertures = CircularAnnulus(positions, r_in=radii[1], r_out=radii[2])

    #Compute raw fluxes
    rawflux_table = aperture_photometry(scidata, apertures)
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
    #print "bkg_mean=",bkg_mean
    bkg_sum = bkg_mean * apertures.area()
    #print "bkg_sum=",bkg_sum
    final_sum = phot_table['aperture_sum_raw'] - bkg_sum
    #print "final_sum=",final_sum
    rate=final_sum/header['EXPTIME']
    phot_table['net_count_rate'] = rate
    #print "Raw=",rawflux_table
    #print 'Bkg=',bkgflux_table
    phot_table['Names']=IDs
    phot_table['Filter']=header['Filter']
    phot_table['Date-Obs']=header['MIDPOINT']
    phot_table.remove_column('id_bkg')
    phot_table.remove_column('xcenter_bkg')
    phot_table.remove_column('ycenter_bkg')
    #print phot_table
    Filter=Meta.FilterParameters(header['FILTER'])
    WVCenter=Filter.CenterWV###Testing Area
    
    pl.figure(figsize=(6,4), dpi=150, facecolor="white")
    pl.imshow(scidata)
    #pl.xlim(0, 170)
    #pl.ylim(130, 250)
    ap_patches = apertures.plot(color='white', lw=0.5,
                           label='Photometry aperture')
    ann_patches = annulus_apertures.plot(color='red', lw=0.5,
                                    label='Background annulus')
    #labels = (ap_patches[0], ann_patches[0])
    #pl.legend(font=10)
    
    return rate,WVCenter,phot_table

def uniform_lat_grid(Latitude,Signal,Fine=False):
    """
    Takes an existing spectrum on a non-standard or even irregular
    wavelength grid and performs linear interpolation to place the data
    on one of four uniform grids:
        1) -100 to 1100nm with 0.5nm bins (includes zeroeth order)
        2) -100 to 1100nm with 0.1nm bins
        3) 115 to 1062.5nm with 0.5nm bins (Pickles standard range)
        4) 115 to 1062.5nm with 0.1nm bins 
    """
    import numpy as np
    from scipy import interpolate

    if Fine:    #Set grid interval
        dlat=1.0
    else:
        dlat=2.0

    LatGrid=np.arange(-90.,90.1,dlat,dtype=float)
    #print Wavelength.size,Signal.size    
    Interp=interpolate.interp1d(Latitude,Signal,kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    SignalonGrid=Interp(LatGrid)

    return LatGrid,SignalonGrid