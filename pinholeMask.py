#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 08:42:17 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import matplotlib.pyplot as plt
#from cnPipeline import *
from helperFunctions import *
#from datetime import datetime
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from numpy.linalg import linalg
from scipy.ndimage import rotate
from scipy.signal import find_peaks_cwt,find_peaks,peak_widths,peak_prominences

def gauss(x,a,x0,sigma,c):
    return c+a*exp(-(x-x0)**2/(2*sigma**2))

def cnFindSpatialPeaks(im,xLoc,peakWidth,avereageFwhm=5):
  
  
  # try to estimate a maximum height from the data
  #TODO: use bad pixel mask to refine this
  if xLoc < 1024:
    maxHeight = np.max(np.mean(im[:,:950],axis=1))
  else:
    maxHeight = np.max(np.mean(im[:,1090:],axis=1))
  
  profile=np.median(im[:,xLoc-avereageFwhm:xLoc+avereageFwhm], axis=1)
  peaks = find_peaks(profile,width=peakWidth, prominence=[50,maxHeight*3],height=[0,maxHeight*3])
  ind = peaks[0] 
  widths = peaks[1]['widths']
  prominences = peaks[1]['prominences']
  heights = peaks[1]['peak_heights']
  return ind, widths, prominences, heights, profile
# total of 42 groups with an extra center pinhole
spacing = 14
center = np.arange(36,2048,47)
left = np.concatenate((center[0:21],center[22:]))+ spacing
right = np.concatenate((center[0:21],center[22:]))- spacing
y = np.arange(2048)

siSpectrum = np.load("data/spectra/modulated-8-SiIX.npy")
tharSpectrum = np.load("data/spectra/SiIX-ThAr-spectrum.npy")
spatialFwhmPinhole = 3 # pixels
spectral = np.concatenate((tharSpectrum,tharSpectrum))

spectral = np.tile(spectral,(2048,1))+0.01

allY = np.concatenate((left,center,right))
profilePinhole = np.zeros(2048)
profileDiskSlit = np.zeros(2048)
profileDiskSlit[512:1536] = profileDiskSlit[512:1536]+1
kernel = gauss(np.arange(100)-50,1,0,20,0)
spatialDisk = np.convolve(profileDiskSlit,kernel,mode="same")/sum(kernel)

for ii in range(len(allY)):
  profilePinhole = profilePinhole + gauss(y,0.1,allY[ii],spatialFwhmPinhole/2.35,0)
imPinhole = np.tile(profilePinhole,(2048,1)).T * spectral
imPinhole = rotate(imPinhole,0)
di = int((imPinhole.shape[0]-2048)/2)
imPinhole = imPinhole[di:di+2048,di:di+2048]
fig, ax=plt.subplots()
ax.imshow(imPinhole,vmin=0,vmax=0.03)


fig, ax=plt.subplots()
ax.plot(tharSpectrum)

fig, ax=plt.subplots()
ax.plot(profilePinhole)

# some real data
#real = fits.open("23marfocus_-40000.fits")[0].data.astype(np.float32)
#ind, widths, prominences, heights, profile= cnFindSpatialPeaks(-1*real,512,[1,10],
#                                                              avereageFwhm=5)
#fig, (ax1,ax2)=plt.subplots(2)
#ax1.plot(y,profile,y[ind],profile[ind],'ro')
#ax1.set_ylim([0,np.max(profile)])
#ax2.plot(ind,widths)
#ax2.set_ylim([1,4])

plt.show()