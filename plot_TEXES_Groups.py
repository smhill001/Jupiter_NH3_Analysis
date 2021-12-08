# -*- coding: utf-8 -*-
"""
Created on Sun Dec 05 08:48:16 2021

@author: Steven Hill
"""

def plot_TEXES_Groups(ax):
    import matplotlib.pyplot as pl
    import scipy
    import numpy as np
    import matplotlib.pyplot as pl
    pth="F:/Astronomy/Projects/SAS 2021 Ammonia/GeminiTEXES2017/ZonalResults/"
    latc = scipy.fromfile(file=pth+"zmean_g1_retnh3_lat.txt", dtype=float, count=-1, sep=" ")
    pressure = scipy.fromfile(file=pth+"zmean_g1_retnh3_pressure.txt", dtype=float, count=-1, sep=" ")
    data=np.zeros((7,61))
    print "data=",data.shape
    PL=23
    Start=61*PL
    End=61*(PL+1)
    for i in range(1,7):
        print i
        tmp = scipy.fromfile(file=pth+"zmean_g"+str(i)+"_retnh3_data.txt", dtype=float, count=-1, sep=" ")
        dat=tmp[Start:End]
        print "dat=",dat.shape, dat
        data[i-1,:]=dat
    print Start,End,pressure[PL]
    latg=Centric_to_Graphic(latc)
    scaled_data_mean=np.mean(data,axis=0)*8.0e4
    scaled_data_std=np.std(data,axis=0)*8.0e4
    ax.plot(latg,scaled_data_mean)
    ax.fill_between(latg, scaled_data_mean-scaled_data_std, scaled_data_mean+scaled_data_std,alpha=.2)
    print data[Start:End]

    
def plot_Teifel(ax):
    import matplotlib.pyplot as pl
    import scipy
    import numpy as np
    from numpy import genfromtxt
    Lats=[-24.0,-12.5,0.0,12.5,24.0]
    EWs=np.array([5.92,6.78,6.75,6.35,5.38])*0.1
    #ax.scatter(Lats,EWs)
    import scipy
    pth="F:/Astronomy/Projects/SAS 2021 Ammonia/"
    Teifel = np.array(genfromtxt(pth+"Teifel2018-Fig7.txt", delimiter=','))
    #TEXESGrid=np.zeros((181,2))
    #latgrid,tmpsig=CNRJ.uniform_lat_grid(TEXESRef[:,0],TEXESRef[:,1],Fine=True)
    #TEXESGrid[:,0]=latgrid[:]
    #TEXESGrid[:,1]=tmpsig[:]
    #print data[:]
    ax.scatter(Teifel[4:24,0],Teifel[4:24,1]*0.1)
    #pl.ylim(0.0,1.2)

def Centric_to_Graphic(Latc):
    import numpy as np
    Req=71492.0
    Rp=66845.0
    Latg=Latc
    for i in range(len(Latc)):
        Latg[i]=np.arctan(((Req/Rp)**2)*np.tan(Latc[i]*np.pi/180.))*180.0/np.pi
    print Latg
    return Latg