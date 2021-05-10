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

#### SET UP INITIAL CONFIGURATION OF FILES (STILL NEED TO ADD OCTOBER)
    
ref_path='F:/Astronomy/Projects/SAS 2021 Ammonia/'
map_path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
filenames=['TEXES-CIRS-blk-TEXES.txt','TEXES-CIRS-red-CIRS.txt',
           'Profile of 2020-09-02 to 15 NH3Abs map stack.csv']

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

ST2000XMRef = genfromtxt(map_path+filenames[2], delimiter=',')
ST2000XMGrid=np.zeros((181,2))
latgrid,tmpsig=CNRJ.uniform_lat_grid(ST2000XMRef[:,0],ST2000XMRef[:,1],Fine=True)
ST2000XMGrid[:,0]=latgrid[:]
ST2000XMGrid[:,1]=tmpsig[:]


pl.figure(figsize=(6.0, 4.0), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=-45
x1=45
xtks=19
y0=0.0
y1=40
ytks=9

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
pl.grid(linewidth=0.2)
pl.tick_params(axis='both', which='major', labelsize=8)
pl.ylabel("Ammonia Mole Fraction",fontsize=8,color="black")
pl.xlabel("Latitude (deg)",fontsize=8)

pl.plot(TEXESGrid[:,0],TEXESGrid[:,1],color='k',label='TEXES',linewidth=0.5)
pl.plot(CIRSGrid[:,0],CIRSGrid[:,1],color='r',label='CIRS',linewidth=0.5)
pl.plot(ST2000XMGrid[:,0],1.5*ST2000XMGrid[:,1]-170.,color='g',label='ST2000XM')
#pl.title(date)
pl.legend()
#AX.plot(latgrid,AvgSignal,color='r',label='NH3/HIA')

pl.savefig('F:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/ProfileComparison.png',dpi=320)
