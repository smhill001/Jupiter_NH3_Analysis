# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:49:14 2018

@author: Steven Hill

PURPOSE:    This code creates arrays of maps or patches to compare data in different
            spectral bands. In particular, it was developed for Jovian NH3 analysis.

HISTORY:    This code evolved from PlanetMapTest.py committed on 2/21/2019 at
            https://github.com/smhill001/PlanetMaps/blob/master/PlanetMapTest.py.
            That code provided simple equatorial and polar map views for any planet
            as well as a bit of limited patch mapping.

"""

def AmmoniaMaps(coords='map',cont=True,DateSelection='All'):
    import sys
    drive='f:'
    sys.path.append(drive+'\\Astronomy\Python Play')
    sys.path.append(drive+'\\Astronomy\Python Play\Util')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
    import matplotlib.path as mpath

    import ConfigFiles as CF
    import PlotUtils as PU
    import cartopy.crs as ccrs
    import scipy.ndimage as nd
    import pylab as pl
    import numpy as np
    import exifread
    from PIL import Image
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
        
    #Data and initialization section
    drive="f:"
    path="/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/"
    
    if DateSelection=='All':
        Dates=["20200720UT","20200729UT","20200730UT","20200731UT",
               "20200902UT","20200903UT","20200904UT",
               "20200913UT","20200914UT","20200915UT",
               "20200924UTa","20200924UTb","20200925UT",
               "20201007UTa","20201007UTb","20201008UT","20201009UT",
               "20210622UT","20210708UT"]#,"20210708UTw"]
    else:
        Dates=DateSelection
    #Dates=["20201009UT","20210622UT","20210708UT","20210708UTw"]
    
    Camera=["CMOS","CMOS","CMOS","CMOS",
            "CCD","CCD","CCD","CCD","CCD","CCD",
            "CCD","CCD","CCD","CCD","CCD","CCD","CCD",
            "CMOS","CMOS","CMOS"]
    
    SessionIDs={'20200720UT':['Jupiter-20200720UT-RGB',
                           'Jupiter-20200720UT-RGBmono',
                           'Jupiter-20200720UT-ClrSlp',
                           'Jupiter-20200720UT-NH3Abs',
                           'Jupiter-20200720UT-889CH4',
                           'Jupiter-20200720UT-380NUV'],
             '20200729UT':['Jupiter-20200729UT-RGB',
                           'Jupiter-20200729UT-RGBmono',
                           'Jupiter-20200729UT-ClrSlp',
                           'Jupiter-20200729UT-NH3Abs',
                           'Jupiter-20200729UT-889CH4',
                           'Jupiter-20200729UT-380NUV'],
             '20200730UT':['Jupiter-20200730UT-RGB',
                           'Jupiter-20200730UT-RGBmono',
                           'Jupiter-20200730UT-ClrSlp',
                           'Jupiter-20200730UT-NH3Abs',
                           'Jupiter-20200730UT-889CH4',
                           'Jupiter-20200730UT-380NUV'],
             '20200731UT':['Jupiter-20200731UT-RGB',
                           'Jupiter-20200731UT-RGBmono',
                           'Jupiter-20200731UT-ClrSlp',
                           'Jupiter-20200731UT-NH3Abs',
                           'Jupiter-20200731UT-889CH4',
                           'Jupiter-20200731UT-380NUV'],
             '20200902UT':['Jupiter-20200902UT-RGB',
                           'Jupiter-20200902UT-RGBmono',
                           'Jupiter-20200902UT-ClrSlp',
                           'Jupiter-20200902UT-NH3Abs',
                           'Jupiter-20200902UT-889CH4'],
             '20200903UT':['Jupiter-20200903UT-RGB',
                           'Jupiter-20200903UT-RGBmono',
                           'Jupiter-20200903UT-ClrSlp',
                           'Jupiter-20200903UT-NH3Abs',
                           'Jupiter-20200903UT-889CH4'],
             '20200904UT':['Jupiter-20200904UT-RGB',
                           'Jupiter-20200904UT-RGBmono',
                           'Jupiter-20200904UT-ClrSlp',
                           'Jupiter-20200904UT-NH3Abs',
                           'Jupiter-20200904UT-889CH4'],
             '20200913UT':['Jupiter-20200913UT-RGB',
                           'Jupiter-20200913UT-RGBmono',
                           'Jupiter-20200913UT-ClrSlp',
                           'Jupiter-20200913UT-NH3Abs',
                           'Jupiter-20200913UT-889CH4'],
             '20200914UT':['Jupiter-20200914UT-RGB',
                           'Jupiter-20200914UT-RGBmono',
                           'Jupiter-20200914UT-ClrSlp',
                           'Jupiter-20200914UT-NH3Abs',
                           'Jupiter-20200914UT-889CH4'],
             '20200915UT':['Jupiter-20200915UT-RGB',
                           'Jupiter-20200915UT-RGBmono',
                           'Jupiter-20200915UT-ClrSlp',
                           'Jupiter-20200915UT-NH3Abs656',
                           'Jupiter-20200915UT-NH3Abs658',
                           'Jupiter-20200915UT-889CH4'],
             '20200924UTa':['Jupiter-20200924UT-RGBa',
                           'Jupiter-20200924UT-RGBmonoa',
                           'Jupiter-20200924UT-ClrSlpa',
                           'Jupiter-20200924UT-NH3Abs656a',
                           'Jupiter-20200924UT-NH3Abs672a',
                           'Jupiter-20200924UT-889CH4a'],
             '20200924UTb':['Jupiter-20200924UT-RGBb',
                           'Jupiter-20200924UT-RGBmonob',
                           'Jupiter-20200924UT-ClrSlpb',
                           'Jupiter-20200924UT-NH3Abs656b',
                           'Jupiter-20200924UT-NH3Abs672b',
                           'Jupiter-20200924UT-889CH4b'],
             '20200925UT':['Jupiter-20200925UT-RGB',
                           'Jupiter-20200925UT-RGBmono',
                           'Jupiter-20200925UT-ClrSlp',
                           'Jupiter-20200925UT-NH3Abs656',
                           'Jupiter-20200925UT-NH3Abs672',
                           'Jupiter-20200925UT-889CH4'],
             '20201007UTa':['Jupiter-20201007UT-RGB',
                           'Jupiter-20201007UT-RGBmono',
                           'Jupiter-20201007UT-ClrSlp',
                           'Jupiter-20201007UT-NH3Abs656a',
                           'Jupiter-20201007UT-NH3Abs672a',
                           'Jupiter-20201007UT-889CH4'],
             '20201007UTb':['Jupiter-20201007UT-RGB',
                           'Jupiter-20201007UT-RGBmono',
                           'Jupiter-20201007UT-ClrSlp',
                           'Jupiter-20201007UT-NH3Abs656b',
                           'Jupiter-20201007UT-NH3Abs672b',
                           'Jupiter-20201007UT-889CH4'],
             '20201008UT':['Jupiter-20201008UT-RGB',
                           'Jupiter-20201008UT-RGBmono',
                           'Jupiter-20201008UT-ClrSlp',
                           'Jupiter-20201008UT-NH3Abs656',
                           'Jupiter-20201008UT-NH3Abs672',
                           'Jupiter-20201008UT-889CH4'],
             '20201009UT':['Jupiter-20201009UT-RGB',
                           'Jupiter-20201009UT-RGBmono',
                           'Jupiter-20201009UT-ClrSlp',
                           'Jupiter-20201009UT-NH3Abs656',
                           'Jupiter-20201009UT-NH3Abs672',
                           'Jupiter-20201009UT-889CH4'],
             '20210622UT':['Jupiter-20210622UT-RGB',
                           'Jupiter-20210622UT-RGBmono',
                           'Jupiter-20210622UT-ClrSlp',
                           'Jupiter-20210622UT-NH3Abs656',
                           'Jupiter-20210622UT-889CH4',
                           'Jupiter-20210622UT-380NUV'],
             '20210708UT':['Jupiter-20210708UT-RGB',
                           'Jupiter-20210708UT-RGBmono',
                           'Jupiter-20210708UT-ClrSlp',
                           'Jupiter-20210708UT-NH3AbsAvg',
                           'Jupiter-20210708UT-889CH4',
                           'Jupiter-20210708UT-380NUV'],
             '20210708UTw':['Jupiter-20210708UT-RGBw',
                           'Jupiter-20210708UT-RGBmonow',
                           'Jupiter-20210708UT-ClrSlpw',
                           'Jupiter-20210708UT-NH3AbsAvgw',
                           'Jupiter-20210708UT-889CH4w',
                           'Jupiter-20210708UT-380NUVw']}
    
    PlotTypes=["a) RGB","b) Reflectivity","c) Continuum Slope",
               "d) NH3","e) 889nm","f) 380nm"]
    stack_CCD=np.zeros((90,90))
    stack_CMOS=np.zeros((90,90))
    stack_ALL=np.zeros((90,90))
    MapSetup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    NH3Setup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    
    #Begin looping over observing sessions
    for Date in Dates:
        First=True
        ytk=True
        iSession=1
        print "SessionIDs[Date][4]=",SessionIDs[Date][3]
        NH3Setup.loadplotparams(drive,SessionIDs[Date][3],"Map")
        testNH3=nd.imread(drive+path+NH3Setup.DataFile,flatten=False)
        testNH3_extent=[int(NH3Setup.X0),int(NH3Setup.X1),int(NH3Setup.Y0),int(NH3Setup.Y1)]
        testNH3_patch=np.copy(testNH3[90-int(NH3Setup.Y1):90-int(NH3Setup.Y0),
                           180+int(NH3Setup.X0):180+int(NH3Setup.X1)])
        kernel = Gaussian2DKernel(1)
        NH3_conv = convolve(testNH3_patch, kernel)
        
        #Begin Looping over individual patches or bands
        for SessionID in SessionIDs[Date]:
            print "SessionID,iSession=",SessionID,iSession
            if First:
                #fig=pl.figure(figsize=(8.0, 2.0), dpi=150, facecolor="white")
                fig=pl.figure(figsize=(3.7,6.5), dpi=150, facecolor="white")
                First=False
            
            MapSetup.loadplotparams(drive,SessionID,"Map")
            print "*********MapSetup.DataFile= ",MapSetup.DataFile
            if iSession % 2 == 0:
                ytk=False
            else:
                ytk=True
            print "Date[0:6]=",Date[0:6]
            xtk=False
            if iSession > 3 and Date[0:6]=="202009":
                xtk=True
            if iSession > 4:
                xtk=True
            #MapSetup.Setup_CaratoPy_Map("PC",1,5,iSession,ytk,ptitle=PlotTypes[iSession-1])
            MapSetup.Setup_CaratoPy_Map("PC",3,2,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
            #if MapSetup.ColorPlane=="Grey":
            #    test=nd.imread(drive+path+MapSetup.DataFile,flatten=True)
            #elif MapSetup.ColorPlane<>"Grey":
            test=nd.imread(drive+path+MapSetup.DataFile,flatten=False)
            
            test_crs = ccrs.PlateCarree()
            #test_extent=[-180,180,int(MapSetup.Y0),int(MapSetup.Y1)]
            test_extent=[int(MapSetup.X0),int(MapSetup.X1),int(MapSetup.Y0),int(MapSetup.Y1)]
            print int(MapSetup.X0),int(MapSetup.X1),int(MapSetup.Y0),int(MapSetup.Y1)
            print "BEFORE CONDITIONAL"
            print test.shape
            if len(test.shape)>2:
                print MapSetup.ColorPlane
                if MapSetup.ColorPlane=="RGB":
                    #test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,:])
                    test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),
                                            180+int(MapSetup.X0):180+int(MapSetup.X1),:])
                    print test_patch.shape
                if MapSetup.ColorPlane=="RED":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,0]
                if MapSetup.ColorPlane=="GRN":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,1]
                if MapSetup.ColorPlane=="BLU":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,2]      
            else:
                if MapSetup.ColorPlane=="Grey":
                    test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),
                                            180+int(MapSetup.X0):180+int(MapSetup.X1)])
                
            print "test_extent=",test_extent
            print "test_patch.shape=",test_patch.shape
            if "889CH4" in SessionID or "380NUV"in SessionID or "RGBmono" in SessionID:
                clrtbl='gist_heat'
                #clrtbl='bwr'
            else:
                clrtbl='gist_heat_r'
                #clrtbl='bwr_r'
                
            if coords=='map':
                pl.imshow(test_patch, clrtbl,origin='upper', transform=ccrs.PlateCarree(), 
                          extent=test_extent,vmin=0,vmax=65000)
                #if iSession >1:
                    #pl.contour(test_patch,origin='upper', transform=ccrs.PlateCarree(), extent=test_extent,
                    #           colors='white', alpha=0.5,levels=np.linspace(25000,65000, 4))
                if cont==True:
                    pl.contour(NH3_conv,origin='upper', transform=ccrs.PlateCarree(), extent=testNH3_extent,
                               colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                               linestyles='solid')
            elif coords=='meridian':
                pl.imshow(test_patch, clrtbl,origin='upper', extent=[-45.,45.,-45.,45.])
            if coords=="20200729":
                pl.imshow(test_patch, clrtbl,origin='upper', extent=[-14.,-104.,-45.,45.])

            pl.subplots_adjust(left=0.10, bottom=0.10, right=1.0, top=0.92,
                        wspace=0.0, hspace=0.25)
            #STACKING OF IMAGES IF WE DO IT IN THE CODE>>>
            if "NH3Abs" in SessionID:
                if test_patch.shape==stack_ALL.shape:
                    stack_ALL=stack_ALL+test_patch
                    if Camera[iSession-1] == "CMOS":
                        stack_CMOS=stack_CMOS+test_patch
                    elif Camera[iSession-1] == "CCD":
                        stack_CCD=stack_CCD+test_patch
            iSession=iSession+1
            
        #fig.tight_layout()
        pl.annotate(Date[0:len(Date)-2], xy=(0.005,0.98),xycoords='figure fraction', horizontalalignment='left', 
                    verticalalignment='top',color='b')
       
        pl.savefig(drive+path+Date+"Jupiter-NH3.png",dpi=300)
    
    ALL_stretch=np.array(255.0*stack_ALL/stack_ALL.max())
    im = Image.fromarray(ALL_stretch.astype(int))
    im.save(drive+path+Date+"Jupiter-NH3-ALL-Data.png")

    CMOS_stretch=np.array(255.0*stack_CMOS/stack_CMOS.max())
    im = Image.fromarray(CMOS_stretch.astype(int))
    im.save(drive+path+Date+"Jupiter-NH3-CMOS-Data.png")

    CCD_stretch=np.array(255.0*stack_CCD/stack_CCD.max())
    im = Image.fromarray(CCD_stretch.astype(int))
    im.save(drive+path+Date+"Jupiter-NH3-CCD-Data.png")

    return 0