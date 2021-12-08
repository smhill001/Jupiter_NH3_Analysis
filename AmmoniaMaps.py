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
    #import matplotlib.path as mpath

    #import ConfigFiles as CF
    import PlotUtils as PU
    import scipy.ndimage as nd
    import pylab as pl
    import numpy as np
    from matplotlib.colors import ListedColormap #to hack a white cmap for the RGB image
    #import exifread
    #import png
    #from PIL import Image
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from datetime import datetime
    import ephem
    import EWLibV006 as EWL
    import plot_TEXES_Groups as PTG
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
    #stack_CCD=np.zeros((90,90))
    #stack_CMOS=np.zeros((90,90))
    #stack_ALL=np.zeros((90,90))
    MeridEWArray=np.zeros((90,len(DateSelection)))
    ZoneEWArray=np.zeros((90,len(DateSelection)))
    #Create PlotSetup objects
    # MapSetup to be looped through for each of the six map raster plots
    # NH3Setup specifically for creation of the contour overlay of NH3 
    #   absorption on all the maps.
    NH3Setup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    #Begin looping over observing sessions
    DateCounter=0
    for Date in DateSelection:
        #Initialize looping indices and flags
        #First=True
        #ytk=True
        #iSession=1
        
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
        
        if int(CM2deg)<10:
            CM2str="00"+str(int(CM2deg))
        elif int(CM2deg)<100:
            CM2str="0"+str(int(CM2deg))
        else:
            CM2str=str(int(CM2deg))

        #Make patch and smooth for NH3Abs contour overlay
        NH3Abs_patch=np.copy(NH3Abs_map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
        if CM2deg<45:
            NH3Abs_patch=np.concatenate((np.copy(NH3Abs_map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                                  np.copy(NH3Abs_map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
        if CM2deg>315:
            NH3Abs_patch=np.concatenate((np.copy(NH3Abs_map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                                         np.copy(NH3Abs_map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
        kernel = Gaussian2DKernel(1)
        print kernel
        NH3Abs_conv = convolve(NH3Abs_patch, kernel)
        print "-------------------> NH3_conv.shape",NH3Abs_conv.shape
        
        ############################
        MeridProfile=np.mean(NH3Abs_patch[:,:],axis=1)
        MeridProfile=MeridProfile*0.2/((2.0**16.)-1.0)+0.9 #Scale to original image ratio range of 0.9 to 1.1.
        MeridTransmission=0.961*(np.flip(MeridProfile,axis=0))
        MeridAbsorption=1.0-MeridTransmission
        MeridEW=0.55/0.039*MeridAbsorption
        MeridProfileStd=np.std(NH3Abs_patch[:,:],axis=1)/((2.0**16.)-1.0)*0.2
        Meriderror=np.flip(MeridProfileStd,axis=0)
        MeridEWerror=0.55/0.039*Meriderror
        Lats=np.linspace(-44.5,44.5,90)
        print DateCounter; Date
        MeridEWArray[:,DateCounter]=MeridEW[:]
        
        ZoneProfile=np.mean(NH3Abs_patch[:,:],axis=0)
        ZoneProfile=ZoneProfile*0.2/((2.0**16.)-1.0)+0.9 #Scale to original image ratio range of 0.9 to 1.1.
        ZoneTransmission=0.961*((ZoneProfile))
        ZoneAbsorption=1.0-ZoneTransmission
        ZoneEW=0.55/0.039*ZoneAbsorption
        ZoneProfileStd=np.std(NH3Abs_patch[:,:],axis=0)/((2.0**16.)-1.0)*0.2
        Zoneerror=np.flip(ZoneProfileStd,axis=0)
        ZoneEWerror=0.55/0.039*Zoneerror
        Lons=np.linspace(-44.5,44.5,90)
        ZoneEWArray[:,DateCounter]=ZoneEW[:]
        
        #print Lats.shape,ProfileTest.shape
        figprof,axsprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
        figprof.suptitle("Ammonia Absorption Profile",x=0.5,ha='center',color='k')
        axsprof[0].plot(Lats,MeridEW)
        axsprof[0].fill_between(Lats, MeridEW-MeridEWerror, MeridEW+MeridEWerror,alpha=.2)
        axsprof[0].grid(linewidth=0.2)
        axsprof[0].xlim=[-45.,45.]
        axsprof[0].ylim=[0.,1.5]
        axsprof[0].set_xticks(np.linspace(-45.,45.,7), minor=False)
        axsprof[0].set_yticks(np.linspace(0.,1.5,7), minor=False)
        axsprof[0].tick_params(axis='both', which='major', labelsize=7)
        axsprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)
        axsprof[0].set_ylabel("Equivalent Width (nm)",fontsize=8)
        axsprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
        
        axsprof[1].plot(Lons,ZoneEW)
        axsprof[1].fill_between(Lons, ZoneEW-ZoneEWerror, ZoneEW+ZoneEWerror,alpha=.2)
        axsprof[1].grid(linewidth=0.2)
        axsprof[1].xlim=[-45.,45.]
        axsprof[1].ylim=[0.,1.5]
        axsprof[1].set_xticks(np.linspace(-45.,45.,7), minor=False)
        axsprof[1].set_yticks(np.linspace(0.,1.5,7), minor=False)
        axsprof[1].tick_params(axis='both', which='major', labelsize=7)
        axsprof[1].set_xlabel("Longitude (deg)",fontsize=8)
        axsprof[1].set_ylabel("Equivalent Width (nm)",fontsize=8)

        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Profile.png",dpi=300)
        ############################
        #Begin Looping over individual patches or bands
        fig,axs=pl.subplots(2,3,figsize=(6.5,4.2), dpi=150, facecolor="white",
                            sharey=True,sharex=True)
        fig.suptitle(Date+", CM2="+str(int(CM2deg))+", ASI120MM",x=0.5,ha='center',color='k')
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
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                           90-LatLims[0]],vmin=0,vmax=65635,
                                           aspect="equal")
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation='vertical',cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        
                        if PlotTypes[iPlot] in ["NH3Abs"]:
                            cbar.ax.set_yticklabels(['1.91', '0.55', '-0.80'])  # vertical colorbar
                        else:
                            cbar.ax.set_yticklabels(['0.9', '1.0', '1.1']) 
                        cbar.ax.tick_params(labelsize=7)#if iSession >1:"""
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
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
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                           90-LatLims[0]],vmin=0,vmax=65635,
                                           aspect="equal")
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation='vertical',cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        #ax=axs[i,j]
                        cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'])  # vertical colorbar
                        cbar.ax.tick_params(labelsize=7)#if iSession >1:"""
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
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
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                           90-LatLims[0]],vmin=0,vmax=65635,
                                           aspect="equal")
                        cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                                           orientation="vertical",cmap='gist_heat',
                                           ax=axs[i,j],fraction=0.046, pad=0.04)
                        #ax=axs[i,j]
                        #cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'],color="w")  # vertical colorbar
                        #cbar.ax.tick_params(labelsize=7,color="w")#if iSession >1:
                        cbar.set_ticks([])
                        cbar.outline.set_visible(False)
                        #cbar.remove()
                        if cont==True:
                            axs[i,j].contour(NH3Abs_conv,origin='upper', extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                                       colors=['w','k'], alpha=0.5,levels=[28000.0,36000.0],linewidths=[0.5,0.5],
                                       linestyles='solid')
                    else:
                        print PlotTypes[iPlot],"off"
                        fig.delaxes(axs[i,j])
        DateCounter=DateCounter+1
        pl.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.92,
                    wspace=0.25, hspace=0.1)
            
            
        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Map.png",dpi=300)
        
    AvgMeridEW=np.mean(MeridEWArray[:,:],axis=1)
    StdMeridEW=np.std(MeridEWArray[:,:],axis=1)
    AvgZoneEW=np.mean(ZoneEWArray[:,:],axis=1)
    StdZoneEW=np.std(ZoneEWArray[:,:],axis=1)
    figavgprof,axsavgprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figavgprof.suptitle("Average Ammonia Absorption Profile",x=0.5,ha='center',color='k')
    axsavgprof[0].plot(Lats,AvgMeridEW)
    axsavgprof[0].fill_between(Lats, AvgMeridEW-StdMeridEW, AvgMeridEW+StdMeridEW,alpha=.2)
    axsavgprof[0].grid(linewidth=0.2)
    axsavgprof[0].xlim=[-45.,45.]
    axsavgprof[0].ylim=[0.,1.2]
    axsavgprof[0].set_xticks(np.linspace(-45.,45.,7), minor=False)
    axsavgprof[0].set_yticks(np.linspace(0.,1.2,7), minor=False)
    axsavgprof[0].tick_params(axis='both', which='major', labelsize=7)
    axsavgprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)
    axsavgprof[0].set_ylabel("Equivalent Width (nm)",fontsize=8)
    axsavgprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
    PTG.plot_Teifel(axsavgprof[0])
    PTG.plot_TEXES_Groups(axsavgprof[0])

    axsavgprof[1].plot(Lons,AvgZoneEW)
    axsavgprof[1].fill_between(Lons, AvgZoneEW-StdZoneEW, AvgZoneEW+StdZoneEW,alpha=.2)
    axsavgprof[1].grid(linewidth=0.2)
    axsavgprof[1].xlim=[-45.,45.]
    axsavgprof[1].ylim=[0.,1.2]
    axsavgprof[1].set_xticks(np.linspace(-45.,45.,7), minor=False)
    axsavgprof[1].set_yticks(np.linspace(0.,1.2,7), minor=False)
    axsavgprof[1].tick_params(axis='both', which='major', labelsize=7)
    axsavgprof[1].set_xlabel("Longitude (deg)",fontsize=8)
    axsavgprof[1].set_ylabel("Equivalent Width (nm)",fontsize=8)
    pl.savefig(drive+path+"NH3 Map Plots/Jupiter-NH3"+"_CMII_"+
               CM2str+"-AvgProfile.png",dpi=300)
        

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