# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:41:07 2021
###############################################################################
NAME:       ProfileComparison.py

PURPOSE:    To plot Jovian meridional profiles extracted manually
            from data mapped and projected in WinJUPOS and compared to 
            mole fraction retrievals from TEXES and Cassini CIRS.
            This code is used to generate Figure TBD for the 2021 SAS
            ammonia project paper. At some point, it might be integrated
            with the Jupiter Bands code(s).
            
INPUTS:     Many CSV Files
            
LIBRARIES:  TBD
                    

###############################################################################
@author: Steven Hill
"""

import sys
drive='f:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry')
sys.path.append(drive+'/Astronomy/Python Play/FITSImageStuff')
sys.path.append(drive+'/Astronomy/Projects/SAS 2021 Project/Analysis')

import scipy
import pylab as pl
import Meta_and_Control_Data_Operations as Meta
import SpecPhotPlot as SPP
from astropy.io import fits
import ComputeNetRateJupiter as CNRJ
from astropy.io import ascii
from astropy.table import Table, hstack, vstack
from os import listdir
import numpy as np
from numpy import genfromtxt
import copy
from scipy import interpolate
import scipy.stats as ST
import matplotlib.image as mpimg

#### SET UP INITIAL CONFIGURATION OF FILES (STILL NEED TO ADD OCTOBER)
    
ref_path='F:/Astronomy/Projects/SAS 2021 Ammonia/'
map_path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
filenames=['TEXES-CIRS-blk-TEXES.txt',
           'TEXES-CIRS-red-CIRS.txt',
           #'Profile of 20200915UTJupiter-NH3-ALL-Data.csv',
           #'Profile of 20200915UTJupiter-NH3-CMOS-Data.csv',
           #'Profile of 20200915UTJupiter-NH3-CCD-Data.csv']
           #'Profile of 20201009UTJupiter-NH3-ALL-Data.csv',
           #'Profile of 20201009UTJupiter-NH3-CMOS-Data.csv',
           #'Profile of 20201009UTJupiter-NH3-CCD-Data.csv']
           'Profile of 20210622UTJupiter-NH3-ALL-Data.csv',
           'Profile of 20210622UTJupiter-NH3-CMOS-Data.csv',
           'Profile of 20210622UTJupiter-NH3-CCD-Data.csv']

bkgimg = mpimg.imread(ref_path+'TEXES-mirrored.JPG')
#### SET UP CUMULATIVE CANVAS AND PLOT


#### LOOP OVER DATES AND COMPUTE RATIOS

#### READ DATA FROM EITHER CONTINUUM OR NH3
TEXESRef = genfromtxt(ref_path+filenames[0], delimiter=',')
TEXESGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(TEXESRef[:,0],TEXESRef[:,1],Fine=True)
TEXESGrid[:,0]=latgrid[:]
TEXESGrid[:,1]=tmpsig[:]

CIRSRef = genfromtxt(ref_path+filenames[1], delimiter=',')
CIRSGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(CIRSRef[:,0],CIRSRef[:,1],Fine=True)
CIRSGrid[:,0]=latgrid[:]
CIRSGrid[:,1]=tmpsig[:]

ALLRef = genfromtxt(map_path+filenames[2], delimiter=',')
ALLRef[:,0]=ALLRef[:,0]-44.5
ALLGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(ALLRef[:,0],ALLRef[:,1],Fine=True)
ALLGrid[:,0]=latgrid[:]
ALLGrid[:,1]=tmpsig[:]

CMOSRef = genfromtxt(map_path+filenames[3], delimiter=',')
CMOSRef[:,0]=CMOSRef[:,0]-44.5
CMOSGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(CMOSRef[:,0],CMOSRef[:,1],Fine=True)
CMOSGrid[:,0]=latgrid[:]
CMOSGrid[:,1]=tmpsig[:]

CCDRef = genfromtxt(map_path+filenames[4], delimiter=',')
CCDRef[:,0]=CCDRef[:,0]-44.5
CCDGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(CCDRef[:,0],CCDRef[:,1],Fine=True)
CCDGrid[:,0]=latgrid[:]
CCDGrid[:,1]=tmpsig[:]

pl.figure(figsize=(5.0, 5.0), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=-45
x1=45
xtks=19
y0=0
y1=40
ytks=9

pl.imshow(bkgimg,extent=[-45.,45.,0.,40.],aspect='auto')
# Set x limits
pl.xlim(x0,x1)
# Set x ticks
#pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
pl.grid(linewidth=0.2)
pl.tick_params(axis='both', which='major', labelsize=12)
pl.ylabel("Ammonia Mole Fraction",fontsize=14,color="black")
pl.xlabel("Latitude (deg)",fontsize=14)

pl.plot(TEXESGrid[:,0],TEXESGrid[:,1],color='k',label='TEXES',linewidth=1)
pl.plot(CIRSGrid[:,0],CIRSGrid[:,1],color='r',label='CIRS',linewidth=1)
#pl.plot(ALLGrid[:,0],ALLGrid[:,1]*0.7+7,color='b',label='Jul 20 - Sep 15',linewidth=3)
pl.plot(ALLGrid[:,0],ALLGrid[:,1]*0.7+5,color='b',label='Jul 20 - Sep 15',linewidth=3)
#pl.plot(CMOSGrid[:,0],CMOSGrid[:,1]*0.7+7,color='g',label='Jul 20-31',linewidth=2)
#pl.plot(CCDGrid[:,0],CCDGrid[:,1]*0.65+5,color='m',label='Sep 2-15',linewidth=2)
#pl.title(date)
pl.legend()
#AX.plot(latgrid,AvgSignal,color='r',label='NH3/HIA')

pl.savefig('F:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/ProfileComparison.png',dpi=320)
