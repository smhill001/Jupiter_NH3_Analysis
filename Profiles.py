# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:41:07 2021
###############################################################################
NAME:       Profiles.py

PURPOSE:    To plot and analyze Jovian meridional profiles extracted manually
            from calibrated ST2000XM images via MaximDL. The initial purpose
            is to look at NH3 absorption, and perhaps CH4 absorption. Thus the 
            inputs are the two profile files (in-band and out-of-band) .
            
INPUTS:     Two CSV files, one in-band and one out-of-band
            
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
    
root_path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'
observations={'20200902UT':{'pixshft':0.1,'latctr':41.8,'rp':36.1,
                            'exp647':10.0,'exp656':10.0},
              '20200903UT':{'pixshft':-0.4,'latctr':52.7,'rp':36.3,
                            'exp647':10.0,'exp656':10.0},
              '20200904UT':{'pixshft':0.0,'latctr':47.0,'rp':36.2,
                            'exp647':14.5, 'exp656':17.0},
              '20200913UT':{'pixshft':0.4,'latctr':44.7,'rp':35.3,
                            'exp647':10.0,'exp656':9.0},
              '20200914UT':{'pixshft':0.4,'latctr':45.0,'rp':35.2,
                            'exp647':3.0,'exp656':4.0},
              '20200915UT':{'pixshft':0.2,'latctr':40.9,'rp':35.1,
                            'exp647':7.5,'exp656':8.0},
              '20200924UT':{'pixshft':-0.2,'latctr':41.7,'rp':34.1,
                            'exp647':8.5,'exp656':9.5},
              '20200925UT':{'pixshft':0.7,'latctr':45.1,'rp':34.0,
                            'exp647':10.0,'exp656':8.5}}
dates=['20200902UT','20200903UT','20200904UT','20200913UT','20200914UT',
       '20200915UT','20200924UT','20200925UT']

#### SET UP CUMULATIVE CANVAS AND PLOT
canvas=pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")
AX=pl.subplot(1, 1, 1)
x0=-60
x1=60
xtks=13
y0=0.95
y1=1.05
ytks=11
# Set x limits
AX.set_xlim(x0,x1)
# Set x ticks
AX.set_xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
AX.set_ylim(y0,y1)
AX.set_yticks(np.linspace(y0,y1,ytks, endpoint=True))
# Set y ticks
AX.grid(linewidth=0.2)
AX.tick_params(axis='both', which='major', labelsize=8)
AX.set_ylabel(r"$I647/I656$",fontsize=8,color="black")
AX.set_xlabel(r"$Latitude (deg)$",fontsize=8)
AX.set_title("Cumulative")

#### LOOP OVER DATES AND COMPUTE RATIOS
First=True
for date in dates:
    path=root_path+date+'/'
    filelist=listdir(path)
    FNList=[]
    for fn in filelist:
        if "csv" in fn:
            if "Aligned" in fn:
                FNList.append(fn)
    print len(FNList)
    print FNList

#### READ DATA FROM EITHER CONTINUUM OR NH3
    for FN in FNList:
        if "656HIA" in FN:
            CNTRaw = genfromtxt(path+FN, delimiter=',')
        if "647CNT" in FN:
            NH3Raw = genfromtxt(path+FN, delimiter=',')
        
    pixelshift=observations[date]['pixshft']
   
    sumratio=np.sum(NH3Raw[:,1]/observations[date]['exp647'])/ \
        np.sum(CNTRaw[:,1]/observations[date]['exp656'])
    print "sumratio=",sumratio

    CNTNrm = copy.deepcopy(CNTRaw)
    CNTNrm[:,1] = CNTRaw[:,1]/np.max(CNTRaw[:,1])
    
    Interp=interpolate.interp1d(CNTNrm[:,0]+pixelshift,CNTNrm[:,1],kind='linear', 
                                copy=True,bounds_error=False, 
                                fill_value=np.NaN,axis=0)  
    SignalonGrid=Interp(CNTNrm[:,0])
    CNTNrm[:,1]=SignalonGrid
    
    NH3Nrm = copy.deepcopy(NH3Raw)
    NH3Nrm[:,1] = NH3Raw[:,1]/np.max(NH3Raw[:,1])
    NH3Abs=copy.deepcopy(NH3Raw)
    NH3Abs[:,1]=NH3Nrm[:,1]/CNTNrm[:,1]
    
    latctr=observations[date]['latctr']
    rp=observations[date]['rp']
    rp=36.1
    Lat=np.arcsin(-(NH3Raw[:,0]-latctr)/rp)*180.0/np.pi
    latgrid,tmpsig=CNRJ.uniform_lat_grid(Lat,NH3Abs[:,1])
    
    if First:
        signalarray=np.zeros([tmpsig.size,1])
        signalarray[:,0]=tmpsig
        First=False
    else:
        signalarray=np.insert(signalarray,1,tmpsig,axis=1)

    
    pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")
    
    pl.subplot(1, 1, 1)
    #Plot Layout Configuration
    x0=0
    x1=100
    xtks=21
    y0=0.0
    y1=1.2
    ytks=13
    
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
    pl.ylabel(r"$Normalized Intensity$",fontsize=8,color="black")
    pl.xlabel(r"$Pixel",fontsize=8)
    
    pl.plot(NH3Nrm[:,0],NH3Nrm[:,1],color='b',label='647NH3')
    pl.plot(CNTNrm[:,0],CNTNrm[:,1],color='r',label='656HIA')
    pl.plot(NH3Abs[:,0],NH3Abs[:,1],color='k',label='NH3/HIA')
    pl.title(date)
    
    pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")
    
    pl.subplot(1, 1, 1)
    #Plot Layout Configuration
    x0=-80
    x1=80
    xtks=17
    y0=0.0
    y1=1.2
    ytks=13
    
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
    pl.ylabel(r"$Albedo$",fontsize=8,color="black")
    pl.xlabel(r"$Latitude (deg)$",fontsize=8)
    
    pl.plot(Lat,NH3Nrm[:,1],color='b',label='647NH3')
    pl.plot(Lat,CNTNrm[:,1],color='r',label='656HIA')
    pl.plot(Lat,NH3Abs[:,1],color='k',label='NH3/HIA')
    pl.title(date)
    
    AX.plot(Lat,NH3Abs[:,1],color='k',label='NH3/HIA',linewidth=0.5)
    
AvgSignal=np.nanmean(signalarray,axis=1)
std=np.nanstd(signalarray,axis=1) 
sem=ST.sem(signalarray,axis=1,ddof=0,nan_policy='omit')
    
MeanSpec=np.zeros([latgrid.size,4])
MeanSpec[:,0]=latgrid
MeanSpec[:,1]=AvgSignal
MeanSpec[:,2]=std
MeanSpec[:,3]=sem

AX.plot(latgrid,AvgSignal,color='r',label='NH3/HIA')

pl.savefig('F:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/Profiles.png',dpi=320)
