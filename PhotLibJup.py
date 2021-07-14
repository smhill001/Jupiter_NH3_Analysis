# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 16:36:08 2021

@author: Steven Hill
"""

def Campaign(dates):
    import sys
    drive='f:'
    sys.path.append(drive+'/Astronomy/Python Play')
    sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry')
    sys.path.append(drive+'/Astronomy/Python Play/FITSImageStuff')
    sys.path.append(drive+'/Astronomy/Projects/SAS 2021 Project/Analysis')
    
    import numpy as np
    from datetime import datetime

    root_path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/'
    pathout='/Astronomy/Projects/SAS 2021 Ammonia/Jupiter_NH3_Analysis/'
    
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
    Names=['0_Jupiter','1_Io','2_Europa','3_Ganymede','4_Callisto',
           'Moons Ratio','Moons StdP','95% Conf','Trans647','NH3 Abs','Trans Conf']
    datetimearray=np.empty([len(dates)],dtype=datetime)
    
    return root_path,pathout,observations,Names,datetimearray

def CreatePhotTable(root_path,pathout,observations,dates,Names,datetimearray):
    
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
    from astropy.table import Table, vstack
    from os import listdir
    import numpy as np
    from datetime import datetime
    
    counter=0
    First1=True
    First2=True
    for date in dates:
    
        # Create plotable date array
        datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')
        counter=counter+1
        
        # Set image source path and get relevant files filtering on co-aligned files
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
    
        # Get location data Jupiter and Moons, only including observable moons for the current Session
        Jupiter=observations[date]['Jupiter']
        JIDs=observations[date]['JIDs']
        moons=observations[date]['moons']
        mIDs=observations[date]['mIDs']
        
        ###########################################################################
        # Loop over files of available observations and compute net rates for the current observing 
        # session. Then add the result row for each FITS file to a summary table.
        First=True
    
        for FN in FNList:
            hdulist=fits.open(path+FN)
            header=hdulist[0].header
            scidata=hdulist[0].data
            moonsradii=[12,20,28]
            Jupiterradii=[50,85,100]
            
            # Compute Jupiter and Moons countrates and concatenate (vertical stack)
            # into a single table
            moonsrate,WVCenter,mtable=CNRJ.ComputeNetRateJupiter(scidata,header,mIDs,date,moons,moonsradii)
            Jupiterrate,WVCenter,jtable=CNRJ.ComputeNetRateJupiter(scidata,header,JIDs,date,Jupiter,Jupiterradii)
            outtable = vstack([jtable, mtable])
                   
            # Create and load summary table. Adds a row for each (filter) observation for 
            # a given observing session (date)
            if First:
                sumtable=outtable
                First=False
            else:
                sumtable=vstack([sumtable,outtable])
            if First2:
                AllTable=outtable
                First2=False
            else:
                AllTable=vstack([AllTable,outtable])
                
        #Write Net Counts summary table each Session with one row per Filter
        oFN=FN[0:10]+'UT-'+'Jupiter-Photometry.csv'
        #os.linsep=''
        ###########################################################################
    
    
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
        print CH4_1000_abs
        """
    ###############################################################################
    
    ascii.write(AllTable,pathout+'AllTable.csv',format='basic',overwrite=True,delimiter=',')
    return AllTable

def SummaryTablePlot(AllTable,dates,Names,datetimearray,MeasFilt,RefFilt):
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
    from astropy.table import Table, vstack
    from os import listdir
    import numpy as np
    from datetime import datetime
    import PhotLibJup as PLJ    
    
    First2=True
    
    for date in dates:
        mask_date=np.array(AllTable['SessionID']==date)
        mask647A=np.array(AllTable['Filter']==MeasFilt)
        mask647date=mask_date & mask647A
        print AllTable[mask647date]
        t647ncA=AllTable[mask647date]['net_count_rate']
        t647dtA=AllTable[mask647date]['Date-Obs']
    
        mask656A=np.array(AllTable['Filter']==RefFilt)
        mask656date=mask_date & mask656A
        t656ncA=AllTable[mask656date]['net_count_rate']
        t656dtA=AllTable[mask656date]['Date-Obs']
        print AllTable[mask656date]
    
        targetsA=AllTable[mask656date]['Names']
    
        print t647ncA,t656ncA
        if len(t647ncA)==len(t656ncA):
            r647A=t647ncA/t656ncA
            r647A.name='ratio_647_over_656'
            print r647A
            rmean_moonsA=np.mean(r647A[1:5])
            rstd_moonsA=np.std(r647A[1:5])
            Conf_moonsA=rstd_moonsA/np.sqrt(np.count_nonzero(r647A[1:5]))
            
            trans647A=r647A[0]/rmean_moonsA
            NH3_absA=1.0-r647A[0]/rmean_moonsA
            Trans_ConfA=Conf_moonsA/rmean_moonsA
            
            testA=Table({'Names':np.array(targetsA),date:np.array(r647A)},names=('Names',date))
        
            if First2:
                YY=Table({'Names':Names})
            YY[date]=0.0
            #print 'len=',len(np.array(r647))
            for n in range(0,len(np.array(r647A))):
                indxA=np.where(YY['Names']==testA['Names'][n])
                #print 'n, indx=',Names[n],n,indx      
                YY[date][indxA]=testA[date][n]
                #XX[date][0]=t647dt#datetime.strptime(t647dt[0],'%Y-%m-%dT%H:%M:%S')
        
            YY[date][5]=rmean_moonsA
            YY[date][6]=rstd_moonsA
            YY[date][7]=Conf_moonsA
            YY[date][8]=trans647A
            YY[date][9]=NH3_absA
            YY[date][10]=Trans_ConfA
        
            First2=False
    
    YY['Mean Ratio']=0.0
    YY['StdP Ratio']=0.0
    YY['Conf 95%']=0.0
    
    pl.figure(figsize=(6,4), dpi=150, facecolor="white")
    print len(dates)
    for i in range(0,6):
        tmparr=np.zeros(len(dates))
        for j in range(0,len(dates)):
            print i,j,dates[j]
            print YY
            tmparr[j]=YY[dates[j]][i]
        tmparr[tmparr == 0] = np.nan
        
        mean_ratio=np.nanmean(tmparr)
        std_ratio=np.nanstd(tmparr)
        Conf95=std_ratio/np.sqrt(np.count_nonzero(tmparr))
        YY['Mean Ratio'][i]=mean_ratio
        YY['StdP Ratio'][i]=std_ratio
        YY['Conf 95%'][i]=Conf95
        #print 'n, indx=',Names[n],n,indx      
        #XX['Mean Ratio'][n]=np.mean(XX[n][1:4])
        starttime=datetime(2020,9,1,0,0,0)
        endtime=datetime(2020,10,10,0,0,0)
    
        pl.xlim(starttime,endtime)
        pl.ylim(1.2,1.45)
        if Names[i]=='0_Jupiter' or Names[i]=='Moons Ratio':
            mkrsize=5.0
        else:
            mkrsize=2.0
        pl.plot_date(datetimearray,tmparr,label=Names[i],xdate=True,fmt='o',markersize=mkrsize)
        pl.legend(fontsize=8)
        pl.grid('both', linewidth=0.5)
        
    return YY