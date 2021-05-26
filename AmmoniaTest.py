# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

This code creates two plots, each with two subplots
    Plot 1:     Presents information on Jupiter reflectance and filters
        Subplot 1:  Karkoschka, 1994 disk-integrated albedo plus a spline fit that removes
                    the NH3 absorption
        Subplot 2:  Filter transmissions computed (647, 656, 658, 672) and convovled
                    with reference albedo (and spline fit for 647)

    Plot 2:     Presents information on potential color biases due to differences in 
                regional/feature albedos on Jupiter and Galilean moon colors
        Subplot 1:  Plots regional reflectivities (I/F) for NEB, EZ, and SEB from Dahl, 2021
        Subplot 2:  Plots Galilean moon reflectivites from Clark & McCord, 1980


@author: Steven Hill
"""
import sys
sys.path.append('f:\\Astronomy\Python Play')
sys.path.append('f:\\Astronomy\Python Play\Util')
sys.path.append('f:\\Astronomy\Python Play\SPLibraries')
sys.path.append('f:\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
from PyAstronomy import pyasl
import GeneralSpecUtils as GSU
import ConfigFiles as CF
import SpecPhotLibV005 as SPL
from numpy import genfromtxt

###### Get reference disk-integrated albedo from Karkoschka, 1994

Jupiter_Karkoschka1993 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Jupiter_Karkoschka1993=scipy.reshape(Jupiter_Karkoschka1993,[Jupiter_Karkoschka1993.size/8,8])

Jupiter_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
Jupiter_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
Jupiter_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,3]

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Jupiter_KarkRef1993[:,0],Jupiter_KarkRef1993[:,1],
                                        Extend=False,Fine=False)

JK=np.zeros((WaveGrid.size,2))
JK[:,0]=WaveGrid
JK[:,1]=SignalonGrid

SplineWV = np.array([560.0, 580.0, 600.0, 635.0, 660.0, 675.0, 690., 714.0])
print SplineWV
Temp=SPL.log_spline_smoother(SplineWV,JK)
Continuum_Albedo=np.zeros((WaveGrid.size,2))
Continuum_Albedo[:,0]=WaveGrid
Continuum_Albedo[:,1]=Temp

CH4_KarkRef1993=np.zeros((Jupiter_Karkoschka1993.size/8,2))
CH4_KarkRef1993[:,0]=Jupiter_Karkoschka1993[:,0]
CH4_KarkRef1993[:,1]=Jupiter_Karkoschka1993[:,2]


###### Get reference regional I/F reflectivities from Dahl, 2021?

Dahl_NEB="F:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-NEB.txt"
NEB = genfromtxt(Dahl_NEB, delimiter=',')
Dahl_SEB="F:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-SEB.txt"
SEB = genfromtxt(Dahl_SEB, delimiter=',')
Dahl_EZ="F:/Astronomy/Projects/SAS 2021 Ammonia/Dahl Spectra-EZ.txt"
EZ = genfromtxt(Dahl_EZ, delimiter=',')

###### Plot disk-integrated reference albedo and simulated ammonia-free albedo

pl.figure(figsize=(4.0, 4.0), dpi=150, facecolor="white")

pl.subplot(2,1, 1)
#Plot Layout Configuration
x0=600.
x1=700.
xtks=6
y0=0.0
y1=1.0
ytks=11

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
pl.ylabel("Albedo",fontsize=8,color="black")
pl.xlabel("Wavelength (nm)",fontsize=8)

pl.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum Albedo',linewidth=1,color='b')
pl.plot(JK[:,0],JK[:,1],label='Jupiter Albedo',linewidth=0.5,color='r')
pl.plot(CH4_KarkRef1993[:,0],CH4_KarkRef1993[:,1],label='CH4 Abs. Coef. ',linewidth=1,color='g')
pl.legend(fontsize=7)


###### Retrieve filter transmissions and convovle with disk integrated albedoes

FilterFile="F:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125024213_1D_WVCal.txt"
FilterOPNC = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNC=scipy.reshape(FilterOPNC,[FilterOPNC.size/2,2])

FilterFile="F:/Astronomy/Projects/Stars/Capella/Spectral Data/1D Spectra/Capella20201125025745_1D_WVCal.txt"
Filter672 = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter672=scipy.reshape(Filter672,[Filter672.size/2,2])
Transmission672=GSU.SpectrumMath(Filter672,FilterOPNC,"Divide")

FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122014325_1D_WVCal.txt"
FilterOPNM = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPNM=scipy.reshape(FilterOPNM,[FilterOPNM.size/2,2])

FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122020503_1D_WVCal.txt"
Filter658 = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter658=scipy.reshape(Filter658,[Filter658.size/2,2])
Transmission658=GSU.SpectrumMath(Filter658,FilterOPNM,"Divide")

FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122015240_1D_WVCal.txt"
Filter656 = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter656=scipy.reshape(Filter656,[Filter656.size/2,2])
Transmission656=GSU.SpectrumMath(Filter656,FilterOPNM,"Divide")

FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122020035_1D_WVCal.txt"
Filter647 = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter647=scipy.reshape(Filter647,[Filter647.size/2,2])
Transmission647=GSU.SpectrumMath(Filter647,FilterOPNM,"Divide")

ContinuumProduct672=GSU.SpectrumMath(Transmission672,Continuum_Albedo,"Multiply")
ContinuumProduct658=GSU.SpectrumMath(Transmission658,Continuum_Albedo,"Multiply")
ContinuumProduct656=GSU.SpectrumMath(Transmission656,Continuum_Albedo,"Multiply")
ContinuumProduct647=GSU.SpectrumMath(Transmission647,Continuum_Albedo,"Multiply")

AbsorptionProduct672=GSU.SpectrumMath(Transmission672,JK,"Multiply")
AbsorptionProduct658=GSU.SpectrumMath(Transmission658,JK,"Multiply")
AbsorptionProduct656=GSU.SpectrumMath(Transmission656,JK,"Multiply")
AbsorptionProduct647=GSU.SpectrumMath(Transmission647,JK,"Multiply")

###### Plot filter transmissions convolved with disk-integrated albedos

pl.subplot(2,1,2)
x0=600.
x1=700.
xtks=6
y0=0.0
y1=0.6
ytks=7

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
pl.ylabel("Albedo x Transmission",fontsize=8,color="black")
pl.xlabel("Wavelength (nm)",fontsize=8)
pl.plot(ContinuumProduct672[:,0],ContinuumProduct672[:,1],label='Continuum Albedo',linewidth=1,color='b')
pl.plot(ContinuumProduct658[:,0],ContinuumProduct658[:,1],linewidth=1,color='b')
pl.plot(ContinuumProduct656[:,0],ContinuumProduct656[:,1],linewidth=1,color='b')
pl.plot(ContinuumProduct647[:,0],ContinuumProduct647[:,1],linewidth=1,color='b')

pl.plot(AbsorptionProduct672[:,0],AbsorptionProduct672[:,1],label='Jupiter Albedo',linewidth=0.5,color='r')
pl.plot(AbsorptionProduct658[:,0],AbsorptionProduct658[:,1],linewidth=0.5,color='r')
pl.plot(AbsorptionProduct656[:,0],AbsorptionProduct656[:,1],linewidth=0.5,color='r')
pl.plot(AbsorptionProduct647[:,0],AbsorptionProduct647[:,1],linewidth=0.5,color='r')

pl.legend(fontsize=7)


pl.savefig('F:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/AmmoniaFilter.png',dpi=320)

###### Plot regional and moon reflectivities

pl.figure(figsize=(4.0, 4.0), dpi=150, facecolor="white")

pl.subplot(2,1, 1)
#Plot Layout Configuration
x0=600.
x1=700.
xtks=6
y0=0.0
y1=1.2
ytks=7

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
pl.ylabel("Albedo",fontsize=8,color="black")
pl.xlabel("Wavelength (nm)",fontsize=8)
"""
pl.plot(NEB[:,0],NEB[:,1],label='Dahl - NEB',linewidth=1,color='r')
pl.plot(SEB[:,0],SEB[:,1]*0.820,label='Dahl - SEB',linewidth=1,color='g')
pl.plot(EZ[:,0],EZ[:,1]*0.715,label='Dahl - EZ',linewidth=1,color='b')

pl.plot(JK[:,0],JK[:,1],label='Kark - regrid',linewidth=1,color='0.3')
pl.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum',linewidth=1,color='0.7')
"""

Callisto1980 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/JovianMoons/References/callisto_no_header.txt", dtype=float, count=-1, sep=" ")    
Callisto1980=scipy.reshape(Callisto1980,[Callisto1980.size/3,3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Callisto1980[:,0]*1000.,Callisto1980[:,1],
                                        Extend=False,Fine=False)
CalGrid=np.zeros((WaveGrid.size,2))
CalGrid[:,0]=WaveGrid
CalGrid[:,1]=SignalonGrid

CallistoProduct672=GSU.SpectrumMath(Transmission672,CalGrid,"Multiply")
CallistoProduct658=GSU.SpectrumMath(Transmission658,CalGrid,"Multiply")
CallistoProduct656=GSU.SpectrumMath(Transmission656,CalGrid,"Multiply")
CallistoProduct647=GSU.SpectrumMath(Transmission647,CalGrid,"Multiply")


Ganymede1980 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/JovianMoons/References/ganymede_no_header.txt", dtype=float, count=-1, sep=" ")    
Ganymede1980=scipy.reshape(Ganymede1980,[Ganymede1980.size/3,3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Ganymede1980[:,0]*1000.,Ganymede1980[:,1],
                                        Extend=False,Fine=False)
GanGrid=np.zeros((WaveGrid.size,2))
GanGrid[:,0]=WaveGrid
GanGrid[:,1]=SignalonGrid

GanymedeProduct672=GSU.SpectrumMath(Transmission672,GanGrid,"Multiply")
GanymedeProduct658=GSU.SpectrumMath(Transmission658,GanGrid,"Multiply")
GanymedeProduct656=GSU.SpectrumMath(Transmission656,GanGrid,"Multiply")
GanymedeProduct647=GSU.SpectrumMath(Transmission647,GanGrid,"Multiply")

Europa1980 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/JovianMoons/References/europa_no_header.txt", dtype=float, count=-1, sep=" ")    
Europa1980=scipy.reshape(Europa1980,[Europa1980.size/3,3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Europa1980[:,0]*1000.,Europa1980[:,1],
                                        Extend=False,Fine=False)
EurGrid=np.zeros((WaveGrid.size,2))
EurGrid[:,0]=WaveGrid
EurGrid[:,1]=SignalonGrid

EuropaProduct672=GSU.SpectrumMath(Transmission672,EurGrid,"Multiply")
EuropaProduct658=GSU.SpectrumMath(Transmission658,EurGrid,"Multiply")
EuropaProduct656=GSU.SpectrumMath(Transmission656,EurGrid,"Multiply")
EuropaProduct647=GSU.SpectrumMath(Transmission647,EurGrid,"Multiply")

Io_leading1980 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/JovianMoons/References/io.leading_no_header.txt", dtype=float, count=-1, sep=" ")    
Io_leading1980=scipy.reshape(Io_leading1980,[Io_leading1980.size/3,3])

Io_trailing1980 = scipy.fromfile(file="F:/Astronomy/Projects/Planets/JovianMoons/References/io.trailing_no_header.txt", dtype=float, count=-1, sep=" ")    
Io_trailing1980=scipy.reshape(Io_trailing1980,[Io_trailing1980.size/3,3])

WaveGrid,SignalonGrid=GSU.uniform_wave_grid(Io_trailing1980[:,0]*1000.,Io_trailing1980[:,1],
                                        Extend=False,Fine=False)
Io_Grid=np.zeros((WaveGrid.size,2))
Io_Grid[:,0]=WaveGrid
Io_Grid[:,1]=SignalonGrid

IoProduct672=GSU.SpectrumMath(Transmission672,Io_Grid,"Multiply")
IoProduct658=GSU.SpectrumMath(Transmission658,Io_Grid,"Multiply")
IoProduct656=GSU.SpectrumMath(Transmission656,Io_Grid,"Multiply")
IoProduct647=GSU.SpectrumMath(Transmission647,Io_Grid,"Multiply")

pl.plot(Callisto1980[:,0]*1000.,Callisto1980[:,1],label='Callisto',linewidth=1,color='b')
pl.plot(Ganymede1980[:,0]*1000.,Ganymede1980[:,1],label='Ganymede',linewidth=1,color='g')
pl.plot(Europa1980[:,0]*1000.,Europa1980[:,1],label='Europa',linewidth=1,color='r')
pl.plot(Io_leading1980[:,0]*1000.,Io_leading1980[:,1],label='Io Leading',linewidth=1,linestyle='dashed',color='k')
pl.plot(Io_trailing1980[:,0]*1000.,Io_trailing1980[:,1],label='Io Trailing',linewidth=1,color='k')
pl.legend(fontsize=7)

pl.subplot(2,1,2)
x0=600.
x1=700.
xtks=6
y0=0.0
y1=1.2
ytks=7

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
pl.ylabel("Albedo x Transmission",fontsize=8,color="black")
pl.xlabel("Wavelength (nm)",fontsize=8)

pl.plot(CallistoProduct672[:,0],CallistoProduct672[:,1],label='Callisto',linewidth=0.5,color='b')
pl.plot(CallistoProduct658[:,0],CallistoProduct658[:,1],linewidth=0.5,color='b')
pl.plot(CallistoProduct656[:,0],CallistoProduct656[:,1],linewidth=0.5,color='b')
pl.plot(CallistoProduct647[:,0],CallistoProduct647[:,1],linewidth=0.5,color='b')

pl.plot(GanymedeProduct672[:,0],GanymedeProduct672[:,1],label='Ganymede',linewidth=0.5,color='g')
pl.plot(GanymedeProduct658[:,0],GanymedeProduct658[:,1],linewidth=0.5,color='g')
pl.plot(GanymedeProduct656[:,0],GanymedeProduct656[:,1],linewidth=0.5,color='g')
pl.plot(GanymedeProduct647[:,0],GanymedeProduct647[:,1],linewidth=0.5,color='g')

pl.plot(EuropaProduct672[:,0],EuropaProduct672[:,1],label='Europa',linewidth=0.5,color='r')
pl.plot(EuropaProduct658[:,0],EuropaProduct658[:,1],linewidth=0.5,color='r')
pl.plot(EuropaProduct656[:,0],EuropaProduct656[:,1],linewidth=0.5,color='r')
pl.plot(EuropaProduct647[:,0],EuropaProduct647[:,1],linewidth=0.5,color='r')

pl.plot(IoProduct672[:,0],IoProduct672[:,1],label='Io Trailing',linewidth=0.5,color='k')
pl.plot(IoProduct658[:,0],IoProduct658[:,1],linewidth=0.5,color='k')
pl.plot(IoProduct656[:,0],IoProduct656[:,1],linewidth=0.5,color='k')
pl.plot(IoProduct647[:,0],IoProduct647[:,1],linewidth=0.5,color='k')

pl.legend(fontsize=7)

pl.savefig('F:/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/ColorSlopes.png',dpi=320)

###############################################################################

StartIndex647=np.where(ContinuumProduct647[:,0]==635.0)
EndIndex647=np.where(ContinuumProduct647[:,0]==659.0)
ContimIntegral647=sum(ContinuumProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
AbsorpIntegral647=sum(AbsorptionProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])

print "########### Jupiter 647NH3/647Cont"
print "index=",StartIndex647[0][0],EndIndex647[0][0]
print "Contin, Absorp=",ContimIntegral647,AbsorpIntegral647
print "Ratio, 1-Ratio=",AbsorpIntegral647/ContimIntegral647,1.0-AbsorpIntegral647/ContimIntegral647
print "1/(1-Ratio)=",1.0/(1.0-AbsorpIntegral647/ContimIntegral647)

StartIndex656=np.where(ContinuumProduct656[:,0]==645.0)
EndIndex656=np.where(ContinuumProduct656[:,0]==672.0)
ContimIntegral656=sum(ContinuumProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
AbsorpIntegral656=sum(AbsorptionProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
print
print "########### Jupiter 656NH3/656Cont"
print "index=",StartIndex656[0][0],EndIndex656[0][0]
print "Contin, Absorp=",ContimIntegral656,AbsorpIntegral656
print "Ratio, 1-Ratio=",AbsorpIntegral656/ContimIntegral656,1.0-AbsorpIntegral656/ContimIntegral656
print "1/(1-Ratio)=",1.0/(1.0-AbsorpIntegral656/ContimIntegral656)
print
print "########### Jupiter 647/656 ratios"
print "647/656 Continuum=",ContimIntegral647/ContimIntegral656
print "647/656 with NH3=",AbsorpIntegral647/AbsorpIntegral656
print "(647/656 with NH3)/(647/656 Continuum)=",(AbsorpIntegral647/AbsorpIntegral656)/(ContimIntegral647/ContimIntegral656)

###############################################################################

CalIntegral647=sum(CallistoProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
GanIntegral647=sum(GanymedeProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
EurIntegral647=sum(EuropaProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])
IoIntegral647=sum(IoProduct647[StartIndex647[0][0]:EndIndex647[0][0],1])

CalIntegral656=sum(CallistoProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
GanIntegral656=sum(GanymedeProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
EurIntegral656=sum(EuropaProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
IoIntegral656=sum(IoProduct656[StartIndex656[0][0]:EndIndex656[0][0],1])
print
print "########### Moons 647/656"
print "Callisto=",CalIntegral647/CalIntegral656
print "Ganymede=",GanIntegral647/GanIntegral656
print "Europa=",EurIntegral647/EurIntegral656
print "Io=",IoIntegral647/IoIntegral656
print
print "########### Moons Correction Factor to Jupiter Continuum Estimate"
print "Callisto=",(ContimIntegral647/ContimIntegral656)/(CalIntegral647/CalIntegral656)
print "Ganymede=",(ContimIntegral647/ContimIntegral656)/(GanIntegral647/GanIntegral656)
print "Europa=",(ContimIntegral647/ContimIntegral656)/(EurIntegral647/EurIntegral656)
print "Io=",(ContimIntegral647/ContimIntegral656)/(IoIntegral647/IoIntegral656)
