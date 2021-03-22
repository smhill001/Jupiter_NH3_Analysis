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

    #Create aperture objects    
    apertures = CircularAperture(positions, r=radii[0])
    annulus_apertures = CircularAnnulus(positions, r_in=radii[1], r_out=radii[2])

    #Compute raw fluxes
    rawflux_table = aperture_photometry(scidata, apertures)
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()
    print "bkg_mean=",bkg_mean
    bkg_sum = bkg_mean * apertures.area()
    print "bkg_sum=",bkg_sum
    final_sum = phot_table['aperture_sum_raw'] - bkg_sum
    print "final_sum=",final_sum
    rate=final_sum/header['EXPTIME']
    phot_table['net_count_rate'] = rate
    print "Raw=",rawflux_table
    print 'Bkg=',bkgflux_table
    phot_table['Names']=IDs
    phot_table['Filter']=header['Filter']
    phot_table['Date-Obs']=header['Date-Obs']
    phot_table.remove_column('id_bkg')
    phot_table.remove_column('xcenter_bkg')
    phot_table.remove_column('ycenter_bkg')
    print phot_table
    WVCenter=Filter.CenterWV###Testing Area
    return rate,WVCenter,phot_table