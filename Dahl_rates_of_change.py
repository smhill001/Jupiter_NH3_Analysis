# -*- coding: utf-8 -*-
"""
Created on Sun Nov 07 15:15:54 2021

@author: Steven Hill
"""
import sys
drive='f:'
sys.path.append(drive+'\\Astronomy\Python Play')
sys.path.append(drive+'\\Astronomy\Python Play\Util')
sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')

import matplotlib.pyplot as pl
import ConfigFiles as CF
import numpy as np

path='F:/Astronomy/Projects/SAS 2021 Ammonia/rate_of_change_arrays/rate_of_change_arrays/'

ammonia=CF.readtextfilelines(path+"ammonia")
arr=np.zeros([8,241])
wave=np.linspace(470,950,241)

pl.figure(figsize=(6.0, 6.0), dpi=150, facecolor="white")

pl.subplot(2,1,1)


for i in range(0,8):
    s=ammonia.CfgLines[i]
    fields=s.split(" ")
    for j in range(0,240):
        arr[i,j]=float(fields[j])
    pl.plot(wave,arr[i,:])
    
ammonia_factor=CF.readtextfilelines(path+"ammonia_factor")
arr_factor=np.zeros([8,241])

pl.subplot(2,1,2)

for i in range(0,8):
    s=ammonia_factor.CfgLines[i]
    fields=s.split(" ")
    for j in range(0,240):
        arr_factor[i,j]=float(fields[j])
    pl.plot(wave,arr_factor[i,:])