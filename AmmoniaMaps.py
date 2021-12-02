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
            
EXAMPLE:    AmmoniaMaps(DateSelection=["2021-09-23"])
"""

def AmmoniaMaps(coords='map',cont=True,DateSelection='All',orientation='Landscape',
                patch_from_config_file=False):
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
    import png
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

    """if DateSelection=='All':
        Dates=["2020-07-20","2020-07-29","2020-07-30","2020-07-31",
               "2020-09-02","2020-09-03","2020-09-04",
               "2020-09-13","2020-09-14","2020-09-15",
               "2020-09-24","2020-09-24","2020-09-25",
               "2020-10-07","2020-10-07","2020-10-08","2020-10-09",
               "2021-06-22","2021-07-08","2021-07-19","2021-07-20",
               "2021-09-10","2021-09-13"]
    else:
        Dates=DateSelection
    #Dates=["20201009UT","20210622UT","20210708UT","20210708"]
    """
    
    Camera=["CMOS","CMOS","CMOS","CMOS",
            "CCD","CCD","CCD","CCD","CCD","CCD",
            "CCD","CCD","CCD","CCD","CCD","CCD","CCD",
            "CMOS","CMOS","CMOS","CMOS",
            "CMOS","CMOS"]
    
    CM2=NewGetMaps()     
     
    print "############################"
    print "len(CM2_L360)=",len(CM2)
    #for i in CM2:
    #    print i
    #print "############################"

    #Set up labels for the six kinds of maps
    PlotTypes=["RGB","Reflectivity","ClrSlp",
               "NH3Abs","889CH4","380NUV"]
    
    #Create empty arrays for the aggregation of NH3 absorption patch data, to
    #  be used as input to ProfileComparison.py
    stack_CCD=np.zeros((90,90))
    stack_CMOS=np.zeros((90,90))
    stack_ALL=np.zeros((90,90))
    
    #Create PlotSetup objects
    # MapSetup to be looped through for each of the six map raster plots
    # NH3Setup specifically for creation of the contour overlay of NH3 
    #   absorption on all the maps.
    NH3Setup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    
    #Begin looping over observing sessions
    for Date in DateSelection:
        #Initialize looping indices and flags
        First=True
        ytk=True
        iSession=1
        
        print "Date=",Date
        #print SessionIDs[Date]
        #print "SessionIDs[Date][3]=",SessionIDs[Date][3]
        #dateconverted=Date[0:4]+"-"+Date[4:6]+"-"+Date[6:8]
        #print "dateconverted=====",dateconverted
        
        MapsforDate=[k for k in CM2 if Date in k]
        
        for i in MapsforDate:
            print i
        
        NH3Abs_fn=[k for k in MapsforDate if "NH3Abs" in k]
        #NH3Abs_map=nd.imread(drive+path+NH3Abs_fn[0],flatten=True)*255
        tmp=load_png(drive+path+NH3Abs_fn[0])
        NH3Abs_map=tmp[:,:,0]

        #testNH3_extent=[int(NH3Setup.X0),int(NH3Setup.X1),int(NH3Setup.Y0),int(NH3Setup.Y1)]
        #print "testNH3_extent=",testNH3_extent
        if patch_from_config_file == True:
            LatLims=[90-int(NH3Setup.Y1),90-int(NH3Setup.Y0)]
            LonLims=[360-int(NH3Setup.X0),360-int(NH3Setup.X1)]
        else:
            strdate=NH3Abs_fn[0][0:15]
            print strdate
            dates=[datetime.strptime(strdate,"%Y-%m-%d-%H%M")]
            date_list_datetime,elev_list,airmass_list,CMI_list,CMII_list,\
                Io_vis_list,Europa_vis_list, \
                Ganymede_vis_list,Callisto_vis_list= \
                EWL.JupiterEphemLists(dates,observer)
            print CMII_list[0]*180./np.pi+45,CMII_list[0]-45*180./np.pi
            LatLims=[45,135]
            CM2deg=CMII_list[0]*180./np.pi
            print "CM2deg=",CM2deg
            LonLims=[360-int(CM2deg+45),360-int(CM2deg-45)]
            print "LonLims=",LonLims
        
        #Make patch and smooth for NH3Abs contour overlay
        NH3Abs_patch=np.copy(NH3Abs_map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
        if CM2deg<45:
            NH3Abs_patch=np.concatenate((np.copy(NH3Abs_map[LatLims[0]:LatLims[1],LonLims[0]:360]),
                                  np.copy(NH3Abs_map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
        if CM2deg>315:
            NH3Abs_patch=np.concatenate((np.copy(NH3Abs_map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                                         np.copy(NH3Abs_map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
        kernel = Gaussian2DKernel(1)
        print kernel
        NH3Abs_conv = convolve(NH3Abs_patch, kernel)
        print "-------------------> NH3_conv.shape",NH3Abs_conv.shape


        #Begin Looping over individual patches or bands
        fig,axs=pl.subplots(2,3,figsize=(6.5,4.2), dpi=150, facecolor="white",
                            sharey=True,sharex=True)
        fig.suptitle(Date+", CM2="+str(int(CM2deg))+", ASI120MM",x=0.02,ha='left',color='b')
        for iPlot in range(0,len(PlotTypes)):
            if orientation=='Landscape':
                i=iPlot/3
                j=np.mod(iPlot,3)
                axs[i,j].grid(linewidth=0.2)
                axs[i,j].ylim=[-45.,45.]
                axs[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
                axs[i,j].set_xticks(np.linspace(360,0,25), minor=False)
                axs[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
                axs[i,j].tick_params(axis='both', which='major', labelsize=7)
                axs[i,j].set_title(PlotTypes[iPlot],fontsize=8)
                #axs[i,j].title.set_size(7)
                #axs[i,j].title.titlepad(0.0)
                print "i,j=", i,j
                #ax=pl.subplot(2, 3, iPlot+1)
                print iPlot,PlotTypes[iPlot]
                if PlotTypes[iPlot] in ["889CH4","380NUV","Reflectivity"]:
                    clrtbl='gist_heat'
                    #clrtbl='bwr'
                else:
                    clrtbl='gist_heat_r'
                    #clrtbl='bwr_r'

                if PlotTypes[iPlot] in ["ClrSlp","NH3Abs","889CH4","380NUV"]:
                    fn=[k for k in MapsforDate if str(PlotTypes[iPlot]) in k]
                    if len(fn) > 0:
                        #jmap=nd.imread(drive+path+fn[0],flatten=True)*255
                        tmp=load_png(drive+path+fn[0])
                        jmap=tmp[:,:,0]
                        print jmap.shape, np.max(jmap)
                        #print jmap2.shape, np.max(jmap2)
                        patch=np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
                        print patch.shape,LonLims[0],LonLims[1]
                        if CM2deg<45:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:360]),
                                                  np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
                        if CM2deg>315:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                                                  np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
                        print patch.shape,LonLims[0],LonLims[1]
   
                        axs[i,j].set_adjustable('box-forced')
                        mapshow=axs[i,j].imshow(patch, clrtbl, origin='upper',  
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],
                                           90-LatLims[1]],vmin=0,vmax=65635,
                                           aspect="equal")
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation='vertical',cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        #ax=axs[i,j]
                        cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                        cbar.ax.tick_params(labelsize=7)#if iSession >1:"""
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],
                                       colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                                       linestyles='solid')
                    else:
                        print PlotTypes[iPlot],"off"
                        fig.delaxes(axs[i,j])
                elif PlotTypes[iPlot] in ["Reflectivity"]:
                    fn=[k for k in MapsforDate if "RGB" in k]
                    if len(fn) > 0:
                        jmap=nd.imread(drive+path+fn[0],flatten=True)*255
                        patch=np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
                        if CM2deg<45:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:360]),
                                                  np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
                        if CM2deg>315:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                                                np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    
                        axs[i,j].set_adjustable('box-forced')
                        mapshow=axs[i,j].imshow(patch, clrtbl,origin='upper',  
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],
                                           90-LatLims[1]],vmin=0,vmax=65635,
                                           aspect="equal")
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation='vertical',cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        #ax=axs[i,j]
                        cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                        cbar.ax.tick_params(labelsize=7)#if iSession >1:"""
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],
                                       colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                                       linestyles='solid')
                    else:
                        print PlotTypes[iPlot],"off"
                        fig.delaxes(axs[i,j])
                else:
                    fn=[k for k in MapsforDate if str(PlotTypes[iPlot]) in k]
                    if len(fn) > 0:
                        jmap=nd.imread(drive+path+fn[0],flatten=False)
                        patch=np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
                        if CM2deg<45:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],LonLims[0]:360]),
                                                  np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
                        if CM2deg>315:
                            patch=np.concatenate((np.copy(jmap[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                                                  np.copy(jmap[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    
                        axs[i,j].set_adjustable('box-forced')
                        mapshow=axs[i,j].imshow(patch, clrtbl,origin='upper',  
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],
                                           90-LatLims[1]],vmin=0,vmax=65635,
                                           aspect="equal")
                        
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation="vertical",cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        #ax=axs[i,j]
                        cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                        cbar.ax.tick_params(labelsize=7)#if iSession >1:
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],
                                       colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                                       linestyles='solid')
                    else:
                        print PlotTypes[iPlot],"off"
                        fig.delaxes(axs[i,j])

        pl.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.92,
                    wspace=0.25, hspace=0.1)
        if int(CM2deg)<10:
            CM2str="00"+str(int(CM2deg))
        elif int(CM2deg)<100:
            CM2str="0"+str(int(CM2deg))
        else:
            CM2str=str(int(CM2deg))
            
            
        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"Jupiter-NH3"+"_CMII_"+
                   CM2str+".png",dpi=300)
        """for ax in axs.flat:
                    ax.label_outer()
        if iPlot == 1 or iPlot == 4: ytk=True
                else: ytk=False
                xtk=False
                if iSession > 3: xtk=True
                MapSetup.Setup_SimpleCyl_Map(2,3,iSession,xtk,ytk,ptitle=PlotTypes[iSession-1])
                ax = pl.subplot(3,2,ptitle=PlotTypes(iPlot))
            elif orientation=='Portrait':
                fig=pl.figure(figsize=(3.7,6.5), dpi=150, facecolor="white") #Landscape
                if iPlot % 2 == 0: ytk=False
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
                    pl.contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[0],90-LatLims[1]],
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
                    pl.contour(NH3Abs_conv,origin='upper', extent=[-45.,45.,-45.,45.],
                               colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                               linestyles='solid')
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
    """
    return 0


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

def load_png(file_path):
    """
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    """
    import png
    import numpy as np

    flow_object = png.Reader(filename=file_path)
    flow_direct = flow_object.asDirect()
    flow_data = list(flow_direct[2])
    (w, h) = flow_direct[3]['size']

    flow = np.zeros((h, w, 3), dtype=np.float64)
    for i in range(len(flow_data)):
        flow[i, :, 0] = flow_data[i][0::3]
        flow[i, :, 1] = flow_data[i][1::3]
        flow[i, :, 2] = flow_data[i][2::3]

    #invalid_idx = (flow[:, :, 2] == 0)
    #flow[:, :, 0:2] = (flow[:, :, 0:2] - 2 ** 15) / 64.0
    #flow[invalid_idx, 0] = 0
    #flow[invalid_idx, 1] = 0

    return flow.astype(np.uint16) 