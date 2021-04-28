# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 12:49:14 2018

@author: Steven Hill
"""

def AmmoniaMaps():
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
    
    drive="f:"
    path="/Astronomy/Projects/Planets/Jupiter/Imaging Data/Mapping/"
    
    Dates=["20200720UT","20200729UT","20200730UT","20200731UT"]
    PlotIDs={'20200720UT':['Jupiter-20200720UT-RGB',
                           'Jupiter-20200720UT-ClrSlp',
                           'Jupiter-20200720UT-NH3Abs',
                           'Jupiter-20200720UT-889CH4',
                           'Jupiter-20200720UT-380NUV'],
             '20200729UT':['Jupiter-20200729UT-RGB',
                           'Jupiter-20200729UT-ClrSlp',
                           'Jupiter-20200729UT-NH3Abs',
                           'Jupiter-20200729UT-889CH4',
                           'Jupiter-20200729UT-380NUV'],
             '20200730UT':['Jupiter-20200730UT-RGB',
                           'Jupiter-20200730UT-ClrSlp',
                           'Jupiter-20200730UT-NH3Abs',
                           'Jupiter-20200730UT-889CH4',
                           'Jupiter-20200730UT-380NUV'],
             '20200731UT':['Jupiter-20200731UT-RGB',
                           'Jupiter-20200731UT-ClrSlp',
                           'Jupiter-20200731UT-NH3Abs',
                           'Jupiter-20200731UT-889CH4',
                           'Jupiter-20200731UT-380NUV'],}
    
    PlotTypes=["R(685)GB","685/550","656/647","889","380"]

    MapSetup=PU.PlotSetup("f:/Astronomy/Python Play/PlanetMaps/MapConfig.txt")
    for Date in Dates:
        First=True
        ytk=True
        iplot=1
        for PlotID in PlotIDs[Date]:
            print "PlotID,iplot=",PlotID,iplot
            if First:
                fig=pl.figure(figsize=(8.0, 2.0), dpi=150, facecolor="white")
                First=False
            
            MapSetup.loadplotparams(drive,PlotID,"Map")
            if iplot>1:
                ytk=False
            MapSetup.Setup_CaratoPy_Map("PC",1,5,iplot,ytk,ptitle=PlotTypes[iplot-1])
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
                    test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),180+int(MapSetup.X0):180+int(MapSetup.X1),:])
                    print test_patch.shape
                if MapSetup.ColorPlane=="RED":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,0]
                if MapSetup.ColorPlane=="GRN":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,1]
                if MapSetup.ColorPlane=="BLU":
                    test_patch=test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),:,2]      
            else:
                if MapSetup.ColorPlane=="Grey":
                    test_patch=np.copy(test[90-int(MapSetup.Y1):90-int(MapSetup.Y0),180+int(MapSetup.X0):180+int(MapSetup.X1)])
                
            print "test_extent=",test_extent
            if "889CH4" in PlotID or "380NUV" in PlotID:
                clrtbl='gist_heat'
            else:
                clrtbl='gist_heat_r'
            pl.imshow(test_patch, clrtbl,origin='upper', transform=ccrs.PlateCarree(), extent=test_extent)
            #if (iplot > 1):
                #pl.yticklabels([])
            pl.subplots_adjust(left=0.08, bottom=0.14, right=0.98, top=0.92,
                        wspace=None, hspace=None)
            iplot=iplot+1
            
        pl.annotate(Date[0:len(Date)-2], xy=(0.005,0.98),xycoords='figure fraction', horizontalalignment='left', 
                    verticalalignment='top',color='b')
       
        pl.savefig(drive+path+Date+"Jupiter-NH3.png",dpi=300)
    
    return 0