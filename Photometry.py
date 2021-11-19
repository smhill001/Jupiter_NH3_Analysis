# -*- coding: utf-8 -*-"""###############################################################################NAME:       Photometry.py - Special Jupiter Version!PURPOSE:    To extract photometric information from multiple targets in             a single FITS image file then plot the results and write the data            to file. This is a variation on the BroadBand_Photometry.py            program in the /Python Play/SpectroPhotometry/Photometry directory.            INPUTS:     A list of sessions (observing dates) and two filters, in-band             (absorption) and out-of-band (continuum) are given            LIBRARIES:  This code calls the SpecPhotLibNew.py library. It is an updated            subset of the SpecPhotLibV006.py library that had grown cumbersome.                    UPDATES:            2021-08-10: This module has had much of its content stripped                        and put into separate library modules. It is now                        essentially a driver program.###############################################################################@author: Steven Hill"""import osos.linesep="\r"import sysdrive='f:'sys.path.append(drive+'/Astronomy/Python Play')sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry')sys.path.append(drive+'/Astronomy/Python Play/FITSImageStuff')sys.path.append(drive+'/Astronomy/Projects/SAS 2021 Project/Analysis')  import pylab as plfrom astropy.io import asciiimport PhotLibJup as PLJfrom datetime import datetimeimport numpy as np###############################################################################DateLims={"2020":[datetime(2020,9,1,0,0,0),datetime(2020,10,10,0,0,0)],          "2021":[datetime(2021,8,1,0,0,0),datetime(2021,9,10,0,0,0)]}SessionDates={"647CNT":{"656HIA":['20200902UT','20200903UT','20200904UT','20200913UT',                           '20200914UT','20200915UT','20200924UT','20200925UT',                           '20201007UT','20201008UT','20201009UT',                           '20210812UT','20210817UT','20210830UT','20210905UT',                           '20210906UT'],                 "672SII":['20200924UT','20200925UT','20201007UT','20201008UT',                           '20201009UT'],                 "658NII":['20200915UT'],                 "632OII":['20210812UT','20210817UT','20210830UT','20210905UT',                           '20210906UT']},       "889CH4":{"940NIR":['20200902UT','20200903UT','20200904UT','20200913UT',                           '20200914UT','20200915UT','20200924UT','20200925UT',                           '20201007UT','20201008UT','20201009UT','20210812UT',                           '20210817UT','20210830UT','20210905UT','20210906UT']}}#******************************************************************************# Initialize metadata for the 2020 647/656 campaign#******************************************************************************MeasFilt='647CNT'RefFilt='656HIA'dates647_656=SessionDates[MeasFilt][RefFilt]dates647_656_2020= [x for x in dates647_656 if x[0:4]=='2020']root_path,pathout,observations=PLJ.Campaign(dates647_656_2020)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates647_656_2020)###############################################################################YY647_656_2020=PLJ.SummaryTablePlot(AllTable,dates647_656_2020,MeasFilt,RefFilt)print YY647_656_2020ascii.write(YY647_656_2020,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)################################################################################******************************************************************************# Initialize metadata for the 2020 647/672 campaign#******************************************************************************MeasFilt='647CNT'RefFilt='672SII'dates647_672=SessionDates[MeasFilt][RefFilt]dates647_672_2020= [x for x in dates647_672 if x[0:4]=='2020']root_path,pathout,observations=PLJ.Campaign(dates647_672_2020)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates647_672_2020)###############################################################################YY647_672_2020=PLJ.SummaryTablePlot(AllTable,dates647_672_2020,MeasFilt,RefFilt)print YY647_672_2020ascii.write(YY647_672_2020,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)################################################################################******************************************************************************# Initialize metadata for the 2020 647/658 campaign#******************************************************************************MeasFilt='647CNT'RefFilt='658NII'dates647_658=SessionDates[MeasFilt][RefFilt]dates647_658_2020= [x for x in dates647_658 if x[0:4]=='2020']root_path,pathout,observations=PLJ.Campaign(dates647_658_2020)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates647_658_2020)###############################################################################YY647_658_2020=PLJ.SummaryTablePlot(AllTable,dates647_658_2020,MeasFilt,RefFilt)print YY647_658_2020ascii.write(YY647_658_2020,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)#******************************************************************************#******************************************************************************#******************************************************************************datetimearray=np.empty([len(dates647_656_2020)],dtype=datetime)tmperr=np.zeros(len(dates647_656_2020))tmparr=np.zeros(len(dates647_656_2020))counter=0for date in dates647_656_2020:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    tmparr[counter]=YY647_656_2020[date][8]    tmperr[counter]=YY647_656_2020[date][10]    counter=counter+1      starttime=DateLims["2020"][0]endtime=DateLims["2020"][1]figsum,ax1=pl.subplots(1,1,figsize=(6,4), dpi=150, facecolor="white")mkrsize=5.0plotshow=ax1.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='647/656',color='k')ax1.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='k',elinewidth=1.0)ax1.plot_date([starttime,endtime],0.961*np.ones(2),xdate=True,             linestyle='dashed',markersize=0.0,color='r',linewidth=1.0)ax1.xlim=[starttime,endtime]ax1.ylim=[0.9,1.05]ax1.yticks=[np.arange(0.90, 1.05, step=0.05)]###############################################################################datetimearray=np.empty([len(dates647_672_2020)],dtype=datetime)tmperr=np.zeros(len(dates647_672_2020))tmparr=np.zeros(len(dates647_672_2020))counter=0for date in dates647_672_2020:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    tmparr[counter]=YY647_672_2020[date][8]    tmperr[counter]=YY647_672_2020[date][10]    counter=counter+1ax1.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='647/672',color='b')ax1.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='b',elinewidth=1.0)###############################################################################datetimearray=np.empty([len(dates647_658_2020)],dtype=datetime)tmperr=np.zeros(len(dates647_658_2020))tmparr=np.zeros(len(dates647_658_2020))counter=0for date in dates647_658_2020:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    tmparr[counter]=YY647_658_2020[date][8]    tmperr[counter]=YY647_658_2020[date][10]    counter=counter+1tmparr[tmparr == 0] = np.nanax1.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='647/658',color='g')ax1.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='g',elinewidth=1.0)###############################################################################ax1.legend(fontsize=8,loc='upper left')ax1.grid('both', linewidth=0.5)pl.savefig(pathout+'Jupiter-Photometry_2020.png',dpi=300)###############################################################################################################################################################******************************************************************************# Initializemetadata for the 2021 647/632 campaign#******************************************************************************dates647_632_2021=['20210812UT','20210817UT','20210830UT','20210905UT','20210906UT']MeasFilt='647CNT'RefFilt='632OI'root_path,pathout,observations=PLJ.Campaign(dates647_632_2021)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates647_632_2021)###############################################################################YY647_632_2021=PLJ.SummaryTablePlot(AllTable,dates647_632_2021,MeasFilt,RefFilt)print YY647_632_2021ascii.write(YY647_632_2021,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)################################################################################******************************************************************************# Initialize metadata for the 2021 647/656 campaign#******************************************************************************MeasFilt='647CNT'RefFilt='656HIA'dates647_656_2021=['20210812UT','20210817UT','20210830UT','20210905UT','20210906UT']root_path,pathout,observations=PLJ.Campaign(dates647_656_2021)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates647_656_2021)###############################################################################YY647_656_2021=PLJ.SummaryTablePlot(AllTable,dates647_656_2021,MeasFilt,RefFilt)print YY647_656_2021ascii.write(YY647_656_2021,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)###############################################################################datetimearray=np.empty([len(dates647_656_2021)],dtype=datetime)counter=0for date in dates647_656_2021:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    counter=counter+1tmperr=np.zeros(len(dates647_656_2021))tmparr=np.zeros(len(dates647_656_2021))  for j in range(0,len(dates647_656_2021)):    print 8,j,dates647_656_2021[j]    #print YY    tmparr[j]=YY647_656_2021[dates647_656_2021[j]][8]    tmperr[j]=YY647_656_2021[dates647_656_2021[j]][10]tmparr[tmparr == 0] = np.nanstarttime=datetime(2021,8,1,0,0,0)endtime=datetime(2021,9,10,0,0,0)pl.figure(figsize=(6,4), dpi=150, facecolor="white")pl.xlim(starttime,endtime)pl.ylim(0.9,1.05)pl.yticks(np.arange(0.90, 1.05, step=0.05))mkrsize=5.0pl.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='647/656',color='k')pl.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='k',elinewidth=1.0)pl.plot_date([starttime,endtime],0.961*np.ones(2),xdate=True,             linestyle='dashed',markersize=0.0,color='r',linewidth=1.0)###############################################################################datetimearray=np.empty([len(dates647_632_2021)],dtype=datetime)counter=0for date in dates647_632_2021:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    counter=counter+1tmperr=np.zeros(len(dates647_632_2021))tmparr=np.zeros(len(dates647_632_2021))  for j in range(0,len(dates647_632_2021)):    print 8,j,dates647_632_2021[j]    #print YY    tmparr[j]=YY647_632_2021[dates647_632_2021[j]][8]    tmperr[j]=YY647_632_2021[dates647_632_2021[j]][10]tmparr[tmparr == 0] = np.nanpl.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='647/632',color='m')pl.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='m',elinewidth=1.0)###############################################################################pl.legend(fontsize=8,loc='upper left')pl.grid('both', linewidth=0.5)pl.savefig(pathout+'Jupiter-Photometry_2021.png',dpi=300)###############################################################################################################################################################******************************************************************************# Initialize metadata for the 2020 889/940 campaign#******************************************************************************dates889_940_2020=['20200902UT','20200903UT','20200904UT','20200913UT','20200914UT',       '20200915UT','20200924UT','20200925UT','20201007UT','20201008UT',       '20201009UT']     MeasFilt='889CH4'RefFilt='940NIR'root_path,pathout,observations=PLJ.Campaign(dates889_940_2020)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates889_940_2020)###############################################################################YY889_940_2020=PLJ.SummaryTablePlot(AllTable,dates889_940_2020,MeasFilt,RefFilt)print YY889_940_2020ascii.write(YY889_940_2020,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)##############################################################################################################################################################datetimearray=np.empty([len(dates889_940_2020)],dtype=datetime)counter=0for date in dates889_940_2020:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    counter=counter+1tmperr=np.zeros(len(dates889_940_2020))tmparr=np.zeros(len(dates889_940_2020))  for j in range(0,len(dates889_940_2020)):    print 8,j,dates889_940_2020[j]    #print YY    tmparr[j]=YY889_940_2020[dates889_940_2020[j]][8]    tmperr[j]=YY889_940_2020[dates889_940_2020[j]][10]tmparr[tmparr == 0] = np.nanstarttime=datetime(2020,9,1,0,0,0)endtime=datetime(2020,10,10,0,0,0)pl.figure(figsize=(6,4), dpi=150, facecolor="white")pl.xlim(starttime,endtime)pl.ylim(0.0,0.201)pl.yticks(np.arange(0.0, 0.201, step=0.05))mkrsize=5.0pl.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='889/940',color='k')pl.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='k',elinewidth=1.0)pl.plot_date([starttime,endtime],0.1*np.ones(2),xdate=True,             linestyle='dashed',markersize=0.0,color='r',linewidth=1.0)###############################################################################pl.legend(fontsize=8,loc='upper left')pl.grid('both', linewidth=0.5)pl.savefig(pathout+'Jupiter-Photometry_2020-CH4.png',dpi=300)###############################################################################################################################################################******************************************************************************# Initialize metadata for the 2021 889/940 campaign#******************************************************************************dates889_940_2021=dates647_656_2021=['20210812UT','20210817UT','20210830UT','20210905UT','20210906UT']     MeasFilt='889CH4'RefFilt='940NIR'root_path,pathout,observations=PLJ.Campaign(dates889_940_2021)################################################################################AllTable contains records for image-target-filter combination (composite key)# Fields include:  id_raw, xcenter_raw, ycenter_raw, #                  aperture_sum_raw, aperture_sum_bkg, net_count_rate, #                  Names, Filter, Date-Obs, SessionIDAllTable=PLJ.CreatePhotTable(root_path,pathout,observations,dates889_940_2021)###############################################################################YY889_940_2021=PLJ.SummaryTablePlot(AllTable,dates889_940_2021,MeasFilt,RefFilt)print YY889_940_2021ascii.write(YY889_940_2021,pathout+'Transmission_'+MeasFilt+'_over_'+RefFilt+'.csv',format='csv',            overwrite=True,delimiter=',')#print pathoutpl.savefig(pathout+'Jupiter-Photometry_'+MeasFilt+'_over_'+RefFilt+'.png',dpi=300)##############################################################################################################################################################datetimearray=np.empty([len(dates889_940_2021)],dtype=datetime)counter=0for date in dates889_940_2021:    # Create plotable date array    datetimearray[counter]=datetime.strptime(date,'%Y%m%dUT')    counter=counter+1tmperr=np.zeros(len(dates889_940_2021))tmparr=np.zeros(len(dates889_940_2021))  for j in range(0,len(dates889_940_2021)):    print 8,j,dates889_940_2021[j]    #print YY    tmparr[j]=YY889_940_2021[dates889_940_2021[j]][8]    tmperr[j]=YY889_940_2021[dates889_940_2021[j]][10]tmparr[tmparr == 0] = np.nanstarttime=datetime(2021,8,1,0,0,0)endtime=datetime(2021,9,10,0,0,0)pl.figure(figsize=(6,4), dpi=150, facecolor="white")pl.xlim(starttime,endtime)pl.ylim(0.0,0.201)pl.yticks(np.arange(0.0, 0.201, step=0.05))mkrsize=5.0pl.plot_date(datetimearray,tmparr,xdate=True,fmt='o',             markersize=mkrsize,label='889/940',color='k')pl.errorbar(datetimearray,tmparr,yerr=tmperr,linewidth=0.0,ecolor='k',elinewidth=1.0)pl.plot_date([starttime,endtime],0.1*np.ones(2),xdate=True,             linestyle='dashed',markersize=0.0,color='r',linewidth=1.0)###############################################################################pl.legend(fontsize=8,loc='upper left')pl.grid('both', linewidth=0.5)pl.savefig(pathout+'Jupiter-Photometry_2021-CH4.png',dpi=300)