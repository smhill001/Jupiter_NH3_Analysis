# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:49:14 2018

@author: Steven Hill

PURPOSE:    This code creates arrays of maps or patches to compare data in 
            different spectral bands. It was developed for Jovian NH3 analysis
            and is designed for six specific maps: RGB, panchromatic
            reflectivity, color slope, NH3 absorption, 889CH4, and 380NUV.
            The default is to plot within 45 degress of latitude and
            longitude of the central meridian and equator.

HISTORY:    This code evolved from PlanetMapTest.py committed on 2/21/2019 at
            https://github.com/smhill001/PlanetMaps/blob/master/PlanetMapTest.py.
            That code provided simple equatorial and polar map views for any planet
            as well as a bit of limited patch mapping.
            
VERSION:    1.0 - Produces all CMOS maps for 2020 and 2021 and all CCD maps
            for 2020. 
            - Has limited or unverified capability to produce cumulative 
              patches for meridional analysis.
            - Does plots Cartopy longitudes, not correct Jovian longitudes
              from WinJUPOS
            - Poorly commented and probably has extraneous code
            - Dependency on F:\Astronomy\Python Play\PlanetMapsMapConfig.txt
            
EXAMPLE:    AmmoniaMaps(DateSelection=["20210923UTa"])
"""

def AmmoniaMaps(coords='map',cont=True,DateSelection='All',orientation='Landscape',
                patch_from_config_file=True):
    import sys
    drive='f:'
    sys.path.append(drive+'\\Astronomy\Python Play')
    sys.path.append(drive+'\\Astronomy\Python Play\Util')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
    sys.path.append(drive+'\\Astronomy\Python Play\SPLibraries')
    import matplotlib.path as mpath

    import ConfigFiles as CF
    import PlotUtils as PU
    import scipy.ndimage as nd
    import pylab as pl
    import numpy as np
    import exifread
    from PIL import Image
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from datetime import datetime
    import ephem
    import EWLibV006 as EWL

    #Data and initialization section
    drive="f:"
    path="/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/"
    
    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')

    if DateSelection=='All':
        Dates=["20200720UT","20200729UT","20200730UT","20200731UT",
               "20200902UT","20200903UT","20200904UT",
               "20200913UT","20200914UT","20200915UT",
               "20200924UTa","20200924UTb","20200925UT",
               "20201007UTa","20201007UTb","20201008UT","20201009UT",
               "20210622UT","20210708UT","20210719UT","20210720UT",
               "20210910UTw","20210913UT"]#,"20210708UTw"]
    else:
        Dates=DateSelection
    #Dates=["20201009UT","20210622UT","20210708UT","20210708UTw"]
    
    Camera=["CMOS","CMOS","CMOS","CMOS",
            "CCD","CCD","CCD","CCD","CCD","CCD",
            "CCD","CCD","CCD","CCD","CCD","CCD","CCD",
            "CMOS","CMOS","CMOS","CMOS",
            "CMOS","CMOS"]
    
    SessionIDs=GetSessionIDs()
    CM2=NewGetMaps()     
     
    print "############################"
    print "len(CM2_L360)=",len(CM2)
    for i in CM2:
        print i
    print "############################"

    #Set up labels for the six kinds of maps
    PlotTypes=["a) RGB","b) Reflectivity","c) Continuum Slope",
               "d) NH3","e) 889nm","f) 380nm"]
    
    #Create empty arrays for the aggregation of NH3 absorption patch data, to
    #  be used as input to ProfileComparison.py
    stack_CCD=np.zeros((90,90))
    stack_CMOS=np.zeros((90,90))
    stack_ALL=np.zeros((90,90))
    
    #Create PlotSetup objects
    # MapSetup to be looped through for each of the six map raster plots
    # NH3Setup specifically for creation of the contour overlay of NH3 
    #   absorption on all the maps.
    MapSetup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    NH3Setup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    
    #Begin looping over observing sessions
    for Date in Dates:
        #Initialize looping indices and flags
        First=True
        ytk=True
        iSession=1
        
        print "Date=",Date
        print SessionIDs[Date]
        print "SessionIDs[Date][3]=",SessionIDs[Date][3]
        dateconverted=Date[0:4]+"-"+Date[4:6]+"-"+Date[6:8]
        print "dateconverted=====",dateconverted
        MapsforDate=[k for k in CM2 if dateconverted in k]
        for i in MapsforDate:
            print i
        NH3=[k for k in MapsforDate if "NH3Abs" in k]
        #Create smoothed contour of NH3 absorption
        NH3Setup.loadplotparams(drive,SessionIDs[Date][3],"Map")
        #testNH3=nd.imread(drive+path+NH3Setup.DataFile,flatten=True)*255
        testNH3=nd.imread(drive+path+NH3[0],flatten=True)*255

        testNH3_extent=[int(NH3Setup.X0),int(NH3Setup.X1),int(NH3Setup.Y0),int(NH3Setup.Y1)]
        print "testNH3_extent=",testNH3_extent
        if patch_from_config_file == True:
            LatLims=[90-int(NH3Setup.Y1),90-int(NH3Setup.Y0)]
            LonLims=[360-int(NH3Setup.X0),360-int(NH3Setup.X1)]
        else:
            strdate=NH3Setup.DataFile[0:15]
            print strdate
            dates=[datetime.strptime(strdate,"%Y-%m-%d-%H%M")]
            date_list_datetime,elev_list,airmass_list,CMI_list,CMII_list,\
                Io_vis_list,Europa_vis_list, \
                Ganymede_vis_list,Callisto_vis_list= \
                EWL.JupiterEphemLists(dates,observer)
            print CMII_list[0]*180./np.pi+45,CMII_list[0]-45*180./np.pi
            LatLims=[45,135]
            LonLims=[360-int(CMII_list[0]*180./np.pi+45),360-int(CMII_list[0]*180./np.pi-45)]
        
        testNH3_patch=np.copy(testNH3[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
              
        kernel = Gaussian2DKernel(1)
        print kernel
        NH3_conv = convolve(testNH3_patch, kernel)
        print "-------------------> NH3_conv.shape",NH3_conv.shape

        #Begin Looping over individual patches or bands
        for SessionID in SessionIDs[Date]:
            #print "SessionID,iSession=",SessionID,iSession
            #Set up figure canvas if the first map of the given session
            if First:
                if orientation=='Landscape':
                    fig=pl.figure(figsize=(6.5,3.7), dpi=150, facecolor="white") #Landscape
                elif orientation=='Portrait':
                    fig=pl.figure(figsize=(3.7,6.5), dpi=150, facecolor="white") #Portrait
                First=False
            #Load map parameters from file F:\Astronomy\Python Play\PlanetMaps\MapConfig.txt
            MapSetup.loadplotparams(drive,SessionID,"Map")
            #Set tick marks according to position in map array
            
            if orientation=='Landscape':
                if iSession == 1 or iSession == 4: ytk=True
                else: ytk=False
                xtk=False
                if iSession > 3: xtk=True
                MapSetup.Setup_SimpleCyl_Map(2,3,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
            elif orientation=='Portrait':
                if iSession % 2 == 0: ytk=False
                else: ytk=True
                xtk=False
                if iSession > 3 and Date[0:6]=="202009": xtk=True
                if iSession > 4: xtk=True
                MapSetup.Setup_SimpleCyl_Map(3,2,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
                
            #Portrait
            #MapSetup.Setup_CaratoPy_Map("PC",3,2,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
            #Landscape
            #MapSetup.Setup_CaratoPy_Map("PC",2,3,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
            #if MapSetup.ColorPlane=="Grey":
            #    test=nd.imread(drive+path+MapSetup.DataFile,flatten=True)
            #elif MapSetup.ColorPlane<>"Grey":
            print "drive+path+MapSetup.DataFile=",drive+path+MapSetup.DataFile
            ###################
            #If I make this conditional on Grey vs RGB, then I could skip a
            #  MaximDL processing step!
            ###################
            print "MapSetup.ColorPlane=",MapSetup.ColorPlane
            if MapSetup.ColorPlane == "Grey":
                test=nd.imread(drive+path+MapSetup.DataFile,flatten=True)*255
            else:
                test=nd.imread(drive+path+MapSetup.DataFile,flatten=False)
            
            #test_extent=[-180,180,int(MapSetup.Y0),int(MapSetup.Y1)]
            test_extent=[int(MapSetup.X0),int(MapSetup.X1),int(MapSetup.Y0),int(MapSetup.Y1)]
            print int(MapSetup.X0),int(MapSetup.X1),int(MapSetup.Y0),int(MapSetup.Y1)
            print "BEFORE CONDITIONAL"
            print test.shape, MapSetup.ColorPlane
            if len(test.shape)>2:
                print MapSetup.ColorPlane
                if MapSetup.ColorPlane=="RGB":
                    #test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,:])
                    test_patch=np.copy(test[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1],:])
                    print "#########Lims:",LatLims[0],LatLims[1],LonLims[0],LonLims[1]
                   # test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),
                   #                         360-int(MapSetup.X0):360-int(MapSetup.X1),:])
                    print "#########Lims:",90-int(MapSetup.Y1),90-int(MapSetup.Y0),\
                        360-int(MapSetup.X0),360-int(MapSetup.X1)
                    
                    print test_patch.shape
                if MapSetup.ColorPlane=="RED":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,0]
                if MapSetup.ColorPlane=="GRN":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,1]
                if MapSetup.ColorPlane=="BLU":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,2]      
            else:
                if MapSetup.ColorPlane=="Grey":
                    #test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),
                    #                        360-int(MapSetup.X0):360-int(MapSetup.X1)])
                    test_patch=np.copy(test[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
                
            print "test_extent=",test_extent
            print "test_patch.shape=",test_patch.shape
            if "889CH4" in SessionID or "380NUV"in SessionID or "RGBmono" in SessionID:
                clrtbl='gist_heat'
                #clrtbl='bwr'
            else:
                clrtbl='gist_heat_r'
                #clrtbl='bwr_r'
                
            if coords=='map':
                mapshow=pl.imshow(test_patch, clrtbl,origin='upper',  
                          extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],vmin=0,vmax=65635)

                cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], orientation='vertical',cmap='gist_heat')
                ax=MapSetup.ax
                cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                cbar.ax.tick_params(labelsize=7)#if iSession >1:
                if cont==True:
                    pl.contour(NH3_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],
                               colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                               linestyles='solid')
            elif coords=='meridian':
                mapshow=pl.imshow(test_patch, clrtbl,origin='upper', extent=[-45.,45.,-45.,45.],vmin=0,vmax=65635)
                if iSession > 1:
                    cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], orientation='vertical',cmap='gist_heat')
                    ax=MapSetup.ax
                    if iSession == 3: cbar.ax.set_yticklabels(['0.5', '1.0', '1.5'])  # vertical colorbar
                    if iSession == 4: cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                     
                    cbar.ax.tick_params(labelsize=6)#if iSession >1:
                if cont==True:
                    pl.contour(NH3_conv,origin='upper', extent=[-45.,45.,-45.,45.],
                               colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                               linestyles='solid')
            if coords=="20200729":
                mapshow=pl.imshow(test_patch, clrtbl,origin='upper', extent=[-14.,-104.,-45.,45.])

            pl.subplots_adjust(left=0.05, bottom=0.10, right=0.95, top=0.92,
                        wspace=0.0, hspace=0.25)
            #STACKING OF IMAGES IF WE DO IT IN THE CODE>>>
            if "NH3Abs" in SessionID:
                print "@@@@@@@@@@@@@@@@ Session ID:",SessionID
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
       
        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"Jupiter-NH3.png",dpi=300)
        
    ###########################################################################
    #Stretch and save average patch data for use in ProfileComparison.py
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

def GetSessionIDs():
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
                       'Jupiter-20210708UT-380NUVw'],
         '20210719UT':['Jupiter-20210719UT-RGB',
                       'Jupiter-20210719UT-RGBmono',
                       'Jupiter-20210719UT-ClrSlp',
                       'Jupiter-20210719UT-NH3AbsAvg',
                       'Jupiter-20210719UT-889CH4',
                       'Jupiter-20210719UT-380NUV'],
         '20210720UT':['Jupiter-20210720UT-RGB',
                       'Jupiter-20210720UT-RGBmono',
                       'Jupiter-20210720UT-ClrSlp',
                       'Jupiter-20210720UT-NH3AbsAvg',
                       'Jupiter-20210720UT-889CH4',
                       'Jupiter-20210720UT-380NUV'],
         '20210910UTw':['Jupiter-20210910UT-RGB',
                       'Jupiter-20210910UT-RGBmono',
                       'Jupiter-20210910UT-ClrSlp',
                       'Jupiter-20210910UT-NH3AbsAvg',
                       'Jupiter-20210910UT-889CH4',
                       'Jupiter-20210910UT-380NUV'],
         '20210913UTw':['Jupiter-20210913UT-RGB',
                       'Jupiter-20210913UT-RGBmono',
                       'Jupiter-20210913UT-ClrSlp',
                       'Jupiter-20210913UT-NH3AbsAvg',
                       'Jupiter-20210913UT-889CH4',
                       'Jupiter-20210913UT-380NUV'],
         '20210915UTw':['Jupiter-20210915UT-RGB',
                       'Jupiter-20210915UT-RGBmono',
                       'Jupiter-20210915UT-ClrSlp',
                       'Jupiter-20210915UT-NH3AbsAvg',
                       'Jupiter-20210915UT-889CH4',
                       'Jupiter-20210915UT-380NUV'],
         '20210919UTw':['Jupiter-20210919UT-RGBw',
                       'Jupiter-20210919UT-RGBmonow',
                       'Jupiter-20210919UT-ClrSlpw',
                       'Jupiter-20210919UT-NH3AbsAvgw'],
         '20210920UTw':['Jupiter-20210920UT-RGB',
                       'Jupiter-20210920UT-RGBmono',
                       'Jupiter-20210920UT-ClrSlp',
                       'Jupiter-20210920UT-NH3AbsAvg',
                       'Jupiter-20210920UT-889CH4',
                       'Jupiter-20210920UT-380NUV'],
         '20210923UTa':['Jupiter-20210923UT-RGBa',
                       'Jupiter-20210923UT-RGBmonoa',
                       'Jupiter-20210923UT-ClrSlpa',
                       'Jupiter-20210923UT-NH3AbsAvga',
                       'Jupiter-20210923UT-889CH4a',
                       'Jupiter-20210923UT-380NUVa'],
          '20210923UTb':['Jupiter-20210923UT-RGBb',
                       'Jupiter-20210923UT-RGBmonob',
                       'Jupiter-20210923UT-ClrSlpb',
                       'Jupiter-20210923UT-NH3AbsAvgb',
                       'Jupiter-20210923UT-889CH4b',
                       'Jupiter-20210923UT-380NUVb'],
          '20210926UTw':['Jupiter-20210926UT-RGBw',
                       'Jupiter-20210926UT-RGBmonow',
                       'Jupiter-20210926UT-ClrSlpw',
                       'Jupiter-20210926UT-NH3AbsAvgw',
                       'Jupiter-20210926UT-380NUVw'],
          '20210927UTw':['Jupiter-20210927UT-RGBw',
                       'Jupiter-20210927UT-RGBmonow',
                       'Jupiter-20210927UT-ClrSlpw',
                       'Jupiter-20210927UT-NH3AbsAvgw',
                       'Jupiter-20210927UT-889CH4w',
                       'Jupiter-20210927UT-380NUVw'],
          '20211017UTw':['Jupiter-20211017UT-RGBw',
                       'Jupiter-20211017UT-RGBmonow',
                       'Jupiter-20211017UT-ClrSlpw',
                       'Jupiter-20211017UT-NH3AbsAvgw',
                       'Jupiter-20211017UT-889CH4w'],
          '20211019UTw':['Jupiter-20211019UT-RGBw',
                       'Jupiter-20211019UT-RGBmonow',
                       'Jupiter-20211019UT-ClrSlpw',
                       'Jupiter-20211019UT-NH3AbsAvgw',
                       'Jupiter-20211019UT-889CH4w'],
          '20211022UTw':['Jupiter-20211022UT-RGBw',
                       'Jupiter-20211022UT-RGBmonow',
                       'Jupiter-20211022UT-ClrSlpw',
                       'Jupiter-20211022UT-NH3AbsAvgw',
                       'Jupiter-20211022UT-889CH4w']}
          
    return SessionIDs

def NewGetMaps():
    import sys
    drive='f:'
    sys.path.append(drive+'\\Astronomy\Python Play')
    sys.path.append(drive+'\\Astronomy\Python Play\Util')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')

    import os
    import scipy.ndimage as nd
    import pylab as pl

    path='F:/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/'
    fnlist = os.listdir(path)
    #print fnlist
    pnglist=[k for k in fnlist if '.png' in k]
    #print pnglist
    CM2_L360=[k for k in pnglist if 'CM2_L360_MAP-BARE' in k]
    #print wavelets
    CM2_L360.sort()
    
    filelistCM2=CM2_L360
    
    return filelistCM2