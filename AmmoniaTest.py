# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 12:38:27 2021

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


pl.figure(figsize=(8.0, 4.0), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=600.
x1=700.
xtks=23
y0=0.0
y1=0.6
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
pl.xlabel(r"$Wavelength (nm)$",fontsize=8)

#pl.plot(Jupiter_KarkRef1993[:,0],Jupiter_KarkRef1993[:,1],label='Karkoschka, 1994',linewidth=1,color='0.5')
#pl.plot(JK[:,0],JK[:,1],label='Kark - regrid',linewidth=1,color='r')


FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122020035_1D_WVCal.txt"
Filter647 = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
Filter647=scipy.reshape(Filter647,[Filter647.size/2,2])

FilterFile="F:/Astronomy/Projects/Planets/Mars/Spectral Data/1D Spectra/Mars20201122014325_1D_WVCal.txt"
FilterOPN = scipy.fromfile(file=FilterFile, dtype=float, count=-1, sep='\t')    
FilterOPN=scipy.reshape(FilterOPN,[FilterOPN.size/2,2])

Transmission=GSU.SpectrumMath(Filter647,FilterOPN,"Divide")
#pl.plot(Transmission[:,0],Transmission[:,1],label='647CNT Transmission',linewidth=1,color='b')

SplineWV = np.array([560.0, 580.0, 600.0, 635.0, 660.0, 675.0, 690., 714.0])
print SplineWV
Temp=SPL.log_spline_smoother(SplineWV,JK)
Continuum_Albedo=np.zeros((WaveGrid.size,2))
Continuum_Albedo[:,0]=WaveGrid
Continuum_Albedo[:,1]=Temp

#pl.plot(Continuum_Albedo[:,0],Continuum_Albedo[:,1],label='Continuum',linewidth=1,color='b')

ContinuumProduct=GSU.SpectrumMath(Transmission,Continuum_Albedo,"Multiply")
AbsorptionProduct=GSU.SpectrumMath(Transmission,JK,"Multiply")

pl.plot(ContinuumProduct[:,0],ContinuumProduct[:,1],label='Continuum',linewidth=1,color='b')
pl.plot(AbsorptionProduct[:,0],AbsorptionProduct[:,1],label='Absorption',linewidth=1,color='b')

StartIndex=np.where(ContinuumProduct[:,0]==635.0)
EndIndex=np.where(ContinuumProduct[:,0]==659.0)
print "index=",StartIndex[0][0],EndIndex[0][0]
ContimIntegral=sum(ContinuumProduct[StartIndex[0][0]:EndIndex[0][0],1])
AbsorpIntegral=sum(AbsorptionProduct[StartIndex[0][0]:EndIndex[0][0],1])
print "Contin, Absorp=",ContimIntegral,AbsorpIntegral
print "Ratio, 1-Ratio=",AbsorpIntegral/ContimIntegral, 1.0-AbsorpIntegral/ContimIntegral
print "1/(1-Ratio)=",1.0/(1.0-AbsorpIntegral/ContimIntegral)