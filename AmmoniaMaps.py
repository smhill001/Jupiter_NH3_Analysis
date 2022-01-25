# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:49:14 2018

@author: Steven Hill

PURPOSE:    This code creates arrays of map patches to evaluate NH3 absorption
            compared to visible RGB imagery, visible panchromatic reflectivity,
            889CH4 band reflectivity, 380NUV band reflectivity, and local 
            color slope.

HISTORY:    This code evolved from PlanetMapTest.py committed on 2/21/2019 at
            https://github.com/smhill001/PlanetMaps/blob/master/PlanetMapTest.py.
            That code provided simple equatorial and polar map views for any planet
            as well as a bit of limited patch mapping.
            
VERSION:    2.0 (12/31/2021) - Produces CMOS and CCD maps and given date as long
            as the file naming and data format conventions are followed.
            - Produces zonal and meridional 645nm EW profiles for each single 
              NH3 observation session.
            - Produces (landscape only) map patch arrays with countour overlays of 
              NH3 absorption EW. Maps patche extent is hardcoded to +/-45 deg
              of the equator and +/-45 deg planetographic latitude from the 
              central meridion of the NH3 absorption map. Maps are computed in
              System II longitude.
            - Produces cumulative (average) profiles with standard deviations
              for a campaign of sessions with overlays of data from Teifel and Fletcher.
            - Deprecates need for an external configuration file (mapconfig.txt).
            
            1.0 - Produces all CMOS maps for 2020 and 2021 and all CCD maps
            for 2020. 
            - Has limited or unverified capability to produce cumulative 
              patches for meridional analysis.
            - Does plots Cartopy longitudes, not correct Jovian longitudes
              from WinJUPOS
            - Poorly commented and probably has extraneous code
            - Dependency on F:\Astronomy\Python Play\PlanetMapsMapConfig.txt
            
EXAMPLE:    AmmoniaMaps(DateSelection=["2021-09-23"])

"""

def AmmoniaMaps(coords='map',cont=True,zonecorr=[0,0],DateSelection='All',orientation='Landscape'):
    ###########################################################################
    # Set up environment and paths
    #
    import sys
    drive='f:'
    sys.path.append(drive+'\\Astronomy\Python Play')
    sys.path.append(drive+'\\Astronomy\Python Play\Util')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
    sys.path.append(drive+'\\Astronomy\Python Play\SPLibraries')
    import scipy.ndimage as nd
    import pylab as pl
    import numpy as np
    from astropy.convolution import Gaussian2DKernel
    from astropy.convolution import convolve
    from datetime import datetime
    import ephem
    import EWLibV006 as EWL
    import plot_TEXES_Groups as PTG

    ###########################################################################
    #Initialization section
    drive="f:"
    path="/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/"
    
    observer = ephem.Observer()
    #Location from Google Maps at 483 S Oneida Way, Denver 80224
    observer.lon = ephem.degrees('-104.907985')
    observer.lat = ephem.degrees('39.708200')
    
    CM2=NewGetMaps()     
    print "############################"
    print "len(CM2_L360)=",len(CM2)

    #Set up labels for the six kinds of maps
    PlotTypes=["NH3Abs","889CH4","380NUV","RGB","Reflectivity","ClrSlp"]
    
    #Create empty arrays for the aggregation of NH3 absorption patch data, to
    #  be used as input to ProfileComparison.py
    MeridEWArray=np.zeros((90,len(DateSelection)))
    ZoneEWArray=np.zeros((90,len(DateSelection)))
    ###########################################################################
    #Begin looping over observing sessions
    DateCounter=0
    for Date in DateSelection:      
        print "Date=",Date     
        ##### Compute Jupiter ephemeris
        MapsforDate=[k for k in CM2 if Date in k]              
        NH3Abs_fn=[k for k in MapsforDate if "NH3Abs" in k]        
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
        
        #Left pad CM2 strings 
        if int(CM2deg)<10:
            CM2str="00"+str(int(CM2deg))
        elif int(CM2deg)<100:
            CM2str="0"+str(int(CM2deg))
        else:
            CM2str=str(int(CM2deg))
            
        #######################################################################
        # Create NH3Abs patch and scale to EW(nm) for contour computations
        tmp=load_png(drive+path+NH3Abs_fn[0])
        NH3Abs_map=tmp[:,:,0]
        NH3Abs_patch=make_patch(NH3Abs_map,LatLims,LonLims,CM2deg)
        NH3Abs_patch_EW=(1.0-0.961*(NH3Abs_patch*0.2/((2.0**16.)-1.0)+0.9))*0.55/0.039
        
        # Make smoothed patch for NH3Abs contours
        kernel = Gaussian2DKernel(1)
        print kernel
        #Empirical longitudinal limb darkening correction
        if zonecorr[:]==[0,0]:
            ZoneLimbDarkCorr=np.ones(90)
        else:
            ZoneLimbDarkCorr=np.mean(NH3Abs_patch_EW[45-zonecorr[0]:45-zonecorr[1],:],axis=0) ###To flip or not to flip in long?
        ZoneLimbDarkCorr=ZoneLimbDarkCorr/np.max(ZoneLimbDarkCorr)
        print ZoneLimbDarkCorr.shape,NH3Abs_patch_EW.shape
        NH3Abs_patch_EW_corr=NH3Abs_patch_EW/ZoneLimbDarkCorr[None,:]

        NH3Abs_conv = convolve(NH3Abs_patch_EW_corr, kernel)
        print "-------------------> NH3_conv.shape",NH3Abs_conv.shape
        
        #######################################################################
        # Compute and plot single-session meridional and zonal profiles from 
        #   unadjusted patch. Also store the profiles for creation of average
        #   profile plots across all sessions in a give run.
        MeridEW=np.flip(np.mean(NH3Abs_patch_EW[:,:],axis=1),axis=0)
        MeridEWerror=np.flip(np.std(NH3Abs_patch_EW[:,:],axis=1),axis=0)
        Lats=np.linspace(-44.5,44.5,90)
        #print DateCounter; Date
        MeridEWArray[:,DateCounter]=MeridEW[:]
        
        ZoneEW=np.mean(NH3Abs_patch_EW[:,:],axis=0) 
        ZoneEWerror=np.std(NH3Abs_patch_EW[:,:],axis=0)
        Lons=np.linspace(-44.5,44.5,90)
        ZoneEWArray[:,DateCounter]=ZoneEW[:]
        
        figprof,axsprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
        figprof.suptitle("Ammonia Absorption Profile",x=0.5,ha='center',color='k')

        axsprof[0].plot(Lats,MeridEW)
        axsprof[0].fill_between(Lats, MeridEW-MeridEWerror, MeridEW+MeridEWerror,alpha=.2)
        axsprof[0].set_title(Date+", CM2="+str(int(CM2deg))+", ASI120MM",fontsize=12)
        axsprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)

        axsprof[1].plot(Lons,ZoneEW)
        axsprof[1].fill_between(Lons, ZoneEW-ZoneEWerror, ZoneEW+ZoneEWerror,alpha=.2)
        axsprof[1].set_xlabel("Longitude from Meridian (deg)",fontsize=8)

        for iprof in range(0,2):
            axsprof[iprof].grid(linewidth=0.2)
            axsprof[iprof].xlim=[-45.,45.]
            axsprof[iprof].ylim=[0.,1.5]
            axsprof[iprof].set_xticks(np.linspace(-45.,45.,7), minor=False)
            axsprof[iprof].set_yticks(np.linspace(0.,1.5,7), minor=False)
            axsprof[iprof].tick_params(axis='both', which='major', labelsize=7)
            axsprof[iprof].set_ylabel("Equivalent Width (nm)",fontsize=8)

        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Profile.png",dpi=300)
        
        #######################################################################
        #Begin Looping over individual patches or bands
        fig,axs=pl.subplots(2,3,figsize=(6.0,4.5), dpi=150, facecolor="white",
                            sharey=True,sharex=True)
        fig.suptitle(Date+", CM2="+str(int(CM2deg)),x=0.5,ha='center',color='k')
        for iPlot in range(0,len(PlotTypes)):
            print "PlotTypes[i]=============",PlotTypes[iPlot]
            if orientation=='Landscape':
                # Set up
                i=iPlot/3                           #Plot row
                j=np.mod(iPlot,3)                   #Plot column
                axs[i,j].grid(linewidth=0.2)
                axs[i,j].ylim=[-45.,45.]
                axs[i,j].xlim=[360-LonLims[0],360-LonLims[1]]
                axs[i,j].set_xticks(np.linspace(450,0,31), minor=False)
                xticklabels=np.array(np.mod(np.linspace(450,0,31),360))
                axs[i,j].set_xticklabels(xticklabels.astype(int))
                axs[i,j].set_yticks(np.linspace(-45,45,7), minor=False)
                axs[i,j].tick_params(axis='both', which='major', labelsize=7)
                print "PlotTypes[iPlot]=============",PlotTypes[iPlot]
                
                if PlotTypes[iPlot] in ["NH3Abs"]:
                    axs[i,j].set_title(PlotTypes[iPlot]+" EW(nm)",fontsize=8)
                else:
                    axs[i,j].set_title(PlotTypes[iPlot],fontsize=8)
                    
                if PlotTypes[iPlot] in ["889CH4","380NUV","Reflectivity"]:
                    clrtbl='gist_heat'
                    #clrtbl='bwr'
                else:
                    clrtbl='gist_heat_r'
                    #clrtbl='bwr_r'

                axs[i,j].set_adjustable('box-forced') 
                
                ###############################################################
                # Read and plot maps patches for each data type
                #
                if PlotTypes[iPlot] in ["ClrSlp","889CH4","380NUV"]:
                    fn=[k for k in MapsforDate if str(PlotTypes[iPlot]) in k]
                    if len(fn) > 0:
                        tmp=load_png(drive+path+fn[0])
                        jmap=tmp[:,:,0]
                        patch=make_patch(jmap,LatLims,LonLims,CM2deg)
   
                        mapshow=axs[i,j].imshow(patch, clrtbl, origin='upper',  
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                           90-LatLims[0]],vmin=0,vmax=65635,
                                           aspect="equal")
                    else:
                        fig.delaxes(axs[i,j])
                        
                elif PlotTypes[iPlot] in ["NH3Abs"]:
                    mapshow=axs[i,j].imshow(NH3Abs_patch_EW_corr, "gist_heat", origin='upper',  
                               extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                       90-LatLims[0]],vmin=0,vmax=1.2,
                                       aspect="equal")
                        
                elif PlotTypes[iPlot] in ["Reflectivity","RGB"]:
                    fn=[k for k in MapsforDate if "RGB" in k]
                    if len(fn) > 0:
                        if PlotTypes[iPlot] in ["Reflectivity"]:
                            jmap=nd.imread(drive+path+fn[0],flatten=True)*255
                            patch=make_patch(jmap,LatLims,LonLims,CM2deg)
                        elif PlotTypes[iPlot] in ["RGB"]:
                            jmap=nd.imread(drive+path+fn[0],flatten=False)
                            patch=make_patch(jmap,LatLims,LonLims,CM2deg)
    
                        mapshow=axs[i,j].imshow(patch, clrtbl,origin='upper',  
                                   extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],
                                           90-LatLims[0]],vmin=0,vmax=65635,
                                           aspect="equal")
                    else:
                        fig.delaxes(axs[i,j])
                # Overplot contours and set axis labels
                if cont==True:
                    temp=make_contours(axs[i,j],NH3Abs_conv,LatLims,LonLims)
                if i==1:
                    axs[i,j].set_xlabel("Sys. II Longitude (deg)",fontsize=8)
                if j==0:
                    axs[i,j].set_ylabel("Planetographic Latitude (deg)",fontsize=8)
                    
        # Increment Date counter, adjust subplots, and save map figure            
        DateCounter=DateCounter+1
        pl.subplots_adjust(left=0.10, bottom=0.08, right=0.98, top=0.90,
                    wspace=0.10, hspace=0.10)            
        pl.savefig(drive+path+"NH3 Map Plots/"+Date+"-Jupiter-NH3"+"_CMII_"+
                   CM2str+"-Map.png",dpi=300)
        
    ###########################################################################
    # Compute and plot cumulative profiles from campaign sessions
    #
    AvgMeridEW=np.mean(MeridEWArray[:,:],axis=1)
    StdMeridEW=np.std(MeridEWArray[:,:],axis=1)
    AvgZoneEW=np.mean(ZoneEWArray[:,:],axis=1)
    StdZoneEW=np.std(ZoneEWArray[:,:],axis=1)
    figavgprof,axsavgprof=pl.subplots(2,1,figsize=(6.0,6.0), dpi=150, facecolor="white")
    figavgprof.suptitle("Average Ammonia Absorption Profile",x=0.5,ha='center',color='k')

    # Meridional profile & reference data
    axsavgprof[0].plot(Lats,AvgMeridEW,label='This Work',color='C0')
    axsavgprof[0].fill_between(Lats, AvgMeridEW-StdMeridEW, AvgMeridEW+StdMeridEW,color='C0',alpha=.2)
    axsavgprof[0].set_xlabel("Planetographic Latitude (deg)",fontsize=8)

    PTG.plot_Teifel(axsavgprof[0],clr='C0')
    axsavgprof[0].legend(fontsize=7,loc=2)

    ax2 = axsavgprof[0].twinx()  # instantiate a second axes that shares the same x-axis
    ax2.ticklabel_format(axis='y',style='sci',scilimits=(0,1), color='C2')
    ax2.tick_params(axis='y', which='major', labelsize=7, color='C2')
    PTG.plot_TEXES_Groups(ax2,clr='C2')

    ax2.set_ylabel("NH3 Mole Fraction at 440mb",fontsize=8)
    ax2.yaxis.label.set_color('C2')
    ax2.legend(fontsize=7,loc=1)

    # Zonal Profile
    axsavgprof[1].plot(Lons,AvgZoneEW)
    axsavgprof[1].fill_between(Lons, AvgZoneEW-StdZoneEW, AvgZoneEW+StdZoneEW,alpha=.2)
    axsavgprof[1].set_xlabel("Longitude (deg)",fontsize=8)
    
    # Plot layout details and labeling
    for iavgprof in range(0,2):
        axsavgprof[iavgprof].grid(linewidth=0.2)
        axsavgprof[iavgprof].xlim=[-45.,45.]
        axsavgprof[iavgprof].ylim=[0.,1.2]
        axsavgprof[iavgprof].set_xticks(np.linspace(-45.,45.,7), minor=False)
        axsavgprof[iavgprof].set_yticks(np.linspace(0.,1.2,7), minor=False)
        axsavgprof[iavgprof].tick_params(axis='both', which='major', labelsize=7)
        axsavgprof[iavgprof].set_ylabel("Equivalent Width (nm)",fontsize=8)
        
    pl.savefig(drive+path+"NH3 Map Plots/Jupiter-NH3_CMII_AvgProfile.png",dpi=300)
        
    return 0


def NewGetMaps():
    import sys
    drive='f:'
    sys.path.append(drive+'\\Astronomy\Python Play')
    sys.path.append(drive+'\\Astronomy\Python Play\Util')
    sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')

    import os

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
    Purpose: Properly load a 48-bit PNG file
    Read from KITTI .png file
    Args:
        file_path string: file path(absolute)
    Returns:
        data (numpy.array): data of image in (Height, Width, 3) layout
    
    FROM: https://www.programcreek.com/python/example/98900/png.Reader
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

    return flow.astype(np.uint16) 

def make_patch(Map,LatLims,LonLims,CM2deg):
    """
    Purpose: Make a map patch and handle the case where the data overlap
             the map edges. This is designed for a map with Jovian longitude
             conventions that with the left boundary at 360 ascending from
             the right boundary at 0. In WinJUPOS, the actual map setting
             shows the left boundary at zero, which is of course, also 360.
    """
    import numpy as np
    patch=np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]:LonLims[1]])
    if CM2deg<45:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],LonLims[0]-1:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]-360])),axis=1)
    if CM2deg>315:
        patch=np.concatenate((np.copy(Map[LatLims[0]:LatLims[1],360+LonLims[0]:360]),
                              np.copy(Map[LatLims[0]:LatLims[1],0:LonLims[1]])),axis=1)
    return patch

def make_contours(ax,NH3Abs_conv,LatLims,LonLims):
    """
    PURPOSE: Overlay countours of NH3 absorption data on Jovian maps.
             Specifically designed for equivalent widths with mean values of
             ~0.55nm
    """
    cs=ax.contour(NH3Abs_conv,origin='upper', 
                  extent=[360-LonLims[0],360-LonLims[1],90-LatLims[1],90-LatLims[0]],
                  colors=['w','w','w','w','w'], alpha=0.5,levels=[0.4,0.5,0.6,0.7,0.8],
                  linewidths=[0.5,0.5,1.0,0.5,0.5],
                  linestyles=['dashed','dashed','solid','dashed','dashed'])
    ax.clabel(cs,[0.4,0.5,0.6,0.7,0.8],inline=True,fmt='%1.1f',fontsize=8)

def make_cbar():
    """
    PURPOSE: This is a collection of code that was once inline when I wanted
             to plot colorbars. Since colorbars were deprecated in the plot
             presentation, I removed the code. But I didn't want to lose it 
             entirely, so I put it here
    """
    """cbar = pl.colorbar(mapshow, ticks=[0.0,0.3, 0.6,0.9,1.2], 
                       orientation='vertical',cmap='gist_heat',
                       ax=axs[i,j],fraction=0.046, pad=0.04)
    
    if PlotTypes[iPlot] in ["NH3Abs"]:
        cbar.ax.set_yticklabels(['0.0','0.3','0.6','0.9','1.2'])  # vertical colorbar
    else:
        cbar.ax.set_yticklabels(['0.9', '1.0', '1.1']) 
    cbar.ax.tick_params(labelsize=7)#if iSession >1:"""
    
    """cbar = pl.colorbar(mapshow, ticks=[0, 32767, 65635], 
                       orientation="vertical",cmap='gist_heat',
                       ax=axs[i,j],fraction=0.046, pad=0.04)
    #ax=axs[i,j]
    #cbar.ax.set_yticklabels(['0.9', '1.0', '1.1'],color="w")  # vertical colorbar
    #cbar.ax.tick_params(labelsize=7,color="w")#if iSession >1:
    cbar.set_ticks([])
    cbar.outline.set_visible(False)
    #cbar.remove()"""