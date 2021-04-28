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
from os import listdir
import numpy as np

root_path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'
observations={'20200902UT':{'mIDs':['1_Io','2_Europa'],
                          'moons':[[1062.392,327.706],[1242.608,333.763]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[879.586,329.834]]},
              '20200903UT':{'mIDs':['3_Ganymede','1_Io'],
                          'moons':[[474.450,704.147],[760.099,696.835]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[988.256,694.938]]},
              '20200904UT':{'mIDs':['3_Ganymede','2_Europa','1_Io','4_Callisto'],
                          'moons':[[313.420,482.301],[507.593,481.448],
                                   [1078.132,484.594],[1404.567,497.557]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[855.706,484.514]]},
              '20200913UT':{'mIDs':['4_Callisto','1_Io','3_Ganymede','2_Europa'],
                          'moons':[[643.549,695.872],[1083.464,719.539],
                                   [1162.836,703.617],[1227.990,723.549]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[926.365,714.827]]},
              '20200914UT':{'mIDs':['2_Europa','1_Io','4_Callisto','3_Ganymede'],
                          'moons':[[649.194,322.115],[758.376,312.837],
                                   [944.940,298.654],[1406.596,313.977]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[870.275,317.233]]},
              '20200915UT':{'mIDs':['2_Europa','4_Callisto','3_Ganymede'],
                          'moons':[[511.397,681.378],
                                   [1163.909,670.529],[1211.903,695.972]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[744.417,688.966]]},
              '20200924UT':{'mIDs':['4_Callisto','3_Ganymede','1_Io','2_Europa'],
                          'moons':[[518.228,711.519],[784.957,706.481],
                                   [1015.554,699.698],[1268.622,703.033]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[1115.138,694.867]]},
              '20200925UT':{'mIDs':['4_Callisto','3_Ganymede','2_Europa','1_Io'],
                          'moons':[[403.446,810.659],[685.120,801.019],
                                   [902.340,799.914],[1385.118,797.310]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[1228.693,800.266]]},
              '20201007UT':{'mIDs':['1_Io','2_Europa','3_Ganymede','4_Callisto'],
                          'moons':[[447.602,750.386],[726.343,747.246],
                                   [854.645,763.049],[1202.839,761.727]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[580.846,753.611]]},
              '20201008UT':{'mIDs':['3_Ganymede','2_Europa','4_Callisto'],
                          'moons':[[627.718,735.023],[1055.585,728.432],
                                   [1125.006,736.472]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[795.991,721.911]]},
              '20201009UT':{'mIDs':['3_Ganymede','2_Europa'],
                          'moons':[[615.799,614.875],[864.000,611.821]],
                          'JIDs':['0_Jupiter'],
                          'Jupiter':[[1098.101,607.210]]}}
dates=['20200902UT','20200903UT','20200904UT','20200913UT','20200914UT',
       '20200915UT','20200924UT','20200925UT','20201007UT','20201008UT',
       '20201009UT']
Names=['0_Jupiter','1_Io','2_Europa','3_Ganymede','4_Callisto',
       'Moons Ratio','Moons StdP','95% Conf','Trans647','NH3 Abs','Trans Conf']

First1=True
for date in dates:
    path=root_path+date+'/'
    filelist=listdir(path)
    FNList=[]
    for fn in filelist:
        if "fit" in fn:
            if "Aligned" in fn:
                FNList.append(fn)
        elif "FITS" in fn:
            if "Aligned" in fn:
                FNList.append(fn)
    print len(FNList)
    print FNList

    moons=observations[date]['moons']
    mIDs=observations[date]['mIDs']
    Jupiter=observations[date]['Jupiter']
    JIDs=observations[date]['JIDs']
    First=True
    
    for FN in FNList:
        #print FN
        hdulist=fits.open(path+FN)
        #print 'HI'
        Filter=Meta.FilterParameters(hdulist[0].header['FILTER'])
        scidata=hdulist[0].data
        header=hdulist[0].header
        moonsradii=[12,20,28]
        Jupiterradii=[50,85,100]
        
        moonsrate,WVCenter,mtable=CNRJ.ComputeNetRateJupiter(scidata,header,mIDs,moons,moonsradii,Filter)
        Jupiterrate,WVCenter,jtable=CNRJ.ComputeNetRateJupiter(scidata,header,JIDs,Jupiter,Jupiterradii,Filter)
        outtable = vstack([jtable, mtable])
        
        pathout='/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/'
        JFN=FN[0:17]+'-'+'Jupiter-'+hdulist[0].header['FILTER']+'-Photometry.csv'
        #print JFN
        #ascii.write(outtable,pathout+JFN,format='csv',overwrite=True,delimiter=',')
        if First:
            #print "First"
            sumtable=outtable
            First=False
        else:
            #print "Second"
            sumtable=vstack([sumtable,outtable])
            
    oFN=FN[0:10]+'UT-'+'Jupiter-Photometry.csv'
    ascii.write(sumtable,pathout+oFN,format='csv',overwrite=True,delimiter=',')
    
    mask647=(sumtable['Filter']=='647CNT')
    #print mask647
    t647nc=sumtable[mask647]['net_count_rate']
    mask656=(sumtable['Filter']=='656HIA')
    #print mask656
    t656nc=sumtable[mask656]['net_count_rate']
    targets=sumtable[mask656]['Names']
    #print t647nc
    #print t656nc
    
    r647=t647nc/t656nc
    r647.name='ratio_647_over_656'
    rmean_moons=np.mean(r647[1:5])
    rstd_moons=np.std(r647[1:5])
    Conf_moons=rstd_moons/np.sqrt(np.count_nonzero(r647[1:5]))

    trans647=r647[0]/rmean_moons
    NH3_abs=1.0-r647[0]/rmean_moons
    Trans_Conf=Conf_moons/rmean_moons
    test=Table({'Names':np.array(targets),date:np.array(r647)},names=('Names',date))
    if First1:
        XX=Table({'Names':Names})
    XX[date]=0.0
    #print 'len=',len(np.array(r647))
    for n in range(0,len(np.array(r647))):
        indx=np.where(XX['Names']==test['Names'][n])
        #print 'n, indx=',Names[n],n,indx      
        XX[date][indx]=test[date][n]

    XX[date][5]=rmean_moons
    XX[date][6]=rstd_moons
    XX[date][7]=Conf_moons
    XX[date][8]=trans647
    XX[date][9]=NH3_abs
    XX[date][10]=Trans_Conf
    
      
    print
    #print test 
    print

    First1=False

    """
    mask940=(sumtable['Filter']=='940NIR')
    t940nc=sumtable[mask940]['net_count_rate']
    mask889=(sumtable['Filter']=='889CH4')
    t889nc=sumtable[mask889]['net_count_rate']
    mask1000=(sumtable['Filter']=='1000NIR')
    t1000nc=sumtable[mask1000]['net_count_rate']
    
    r889=t889nc/t940nc
    trans889=r889[0]/np.mean(r889[1:5])
    CH4_889_abs=1.0-r889[0]/np.mean(r889[1:5])
    
    print
    print r889
    print
    print trans889
    print CH4_889_abs

    r1000=t1000nc/t940nc
    trans1000=r1000[0]/np.mean(r1000[1:5])
    CH4_1000_abs=1.0-r1000[0]/np.mean(r1000[1:5])
  
    print
    print r1000
    print
    print trans1000
    print CH4_1000_abs"""

XX['Mean Ratio']=0.0
XX['StdP Ratio']=0.0
XX['Conf 95%']=0.0
for i in range(0,5):
    tmparr=np.zeros(len(dates))
    for j in range(0,len(dates)):
        tmparr[j]=XX[dates[j]][i]
    tmparr[tmparr == 0] = np.nan
    mean_ratio=np.nanmean(tmparr)
    std_ratio=np.nanstd(tmparr)
    Conf95=std_ratio/np.sqrt(np.count_nonzero(tmparr))
    XX['Mean Ratio'][i]=mean_ratio
    XX['StdP Ratio'][i]=std_ratio
    XX['Conf 95%'][i]=Conf95
    #print 'n, indx=',Names[n],n,indx      
    #XX['Mean Ratio'][n]=np.mean(XX[n][1:4])


print XX
ascii.write(XX,pathout+'Transmission.csv',format='csv',overwrite=True,delimiter=',')
