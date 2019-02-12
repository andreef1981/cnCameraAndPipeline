#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 09:47:51 2019

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
from helperFunctions import *
import matplotlib.pyplot as plt
from scipy import exp


def cnGauss2d(x,y,x0,y0,sx,sy):
  mx, my = np.meshgrid(x, y)
  
  return exp(-(((mx-x0)**2/(2*sx**2))+((my-y0)**2/(2*sy**2))))


# context imager
xscale = 20 #pixel/arcsec
yscale = 20 #pixel/arcsec
xsize = 2048
ysize = 2048
width = 10.

x = 890
y = 1024
widths = np.array([24,22,20,18,15,10,16,19,23,25])
#%%
ciIm = np.zeros((len(widths),xsize,ysize))
for i in range(len(widths)):
  ciIm[i,:,:] = cnGauss2d(np.arange(xsize), np.arange(ysize), x,
      y, widths[i],widths[i])
np.save("data/spectra/defocused-10-CI.npy", ciIm)  
 
#%%  
im =ciIm 
fig, ax=plt.subplots(num=2)
ax.imshow((im[0,:,:]))



