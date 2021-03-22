# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 09:41:47 2017
###############################################################################
NAME:       Photometry.py

PURPOSE:    To extract photometric information from multiple targets in 
            a single FITS image file then plot the results and write the data
            to file. This is a variation on the BroadBand_Photometry.py
            program in the /Python Play/SpectroPhotometry/Photometry directory.
            
INPUTS:     A single parameter, "Target", points to a configuration file
            that provides the stellar targets and plotting information. That
            configuration file also points to secondary configuraiton files
            that contain the FITS data file lists for each observation.
            
LIBRARIES:  This code calls the SpecPhotLibNew.py library. It is an updated
            subset of the SpecPhotLibV006.py library that had grown cumbersome.
                    

###############################################################################
@author: Steven Hill
"""
import sys
drive='f:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry')
sys.path.append(drive+'/Astronomy/Python Play/FITSImageStuff')
sys.path.append(drive+'/Astronomy/Projects/SAS 2021 Project/Analysis')

import pylab as pl
import Meta_and_Control_Data_Operations as Meta
import SpecPhotPlot as SPP
from astropy.io import fits
import ComputeNetRateJupiter as CNRJ
from astropy.io import ascii
from astropy.table import Table, hstack, vstack


path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/20200913UT/'
FN='2020-09-13-0327_1-Jupiter-647CNT-sum10s-Aligned.fit'
moons=[[643.549,695.872],[1083.464,719.539],
           [1162.836,703.617],[1227.990,723.549]]
mIDs=['Callisto','Io','Ganymede','Europa']
Jupiter=[[926.365,714.827]]
JIDs=['Jupiter']
hdulist=fits.open(path+FN)
print 'HI'
Filter=Meta.FilterParameters(hdulist[0].header['FILTER'])
scidata=hdulist[0].data
header=hdulist[0].header
moonsradii=[10,20,30]
Jupiterradii=[40,80,100]

moonsrate,WVCenter,mtable=ComputeNetRateJupiter(scidata,header,mIDs,moons,moonsradii)
Jupiterrate,WVCenter,jtable=ComputeNetRateJupiter(scidata,header,JIDs,Jupiter,Jupiterradii)
outtable = vstack([jtable, mtable])

pathout='/Astronomy/Projects/SAS 2021 Project/Analysis/'
JFN=FN[0:17]+'-'+'Jupiter-'+hdulist[0].header['FILTER']+'-Photometry.csv'
print JFN
ascii.write(outtable,pathout+JFN,format='csv',overwrite=True,delimiter=',')

