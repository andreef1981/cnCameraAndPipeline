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

def cnGauss(x,a,x0,sigma,c):
    return c+a*exp(-(x-x0)**2/(2*sigma**2))

def cnFindSpatialPeaks(im,peakWidth,xLoc=[0,950],averageFwhm=5):
  
  # try to estimate a maximum height from the data
  #TODO: use bad pixel mask to refine this
  #TODO: columns for left and right side are hard coded
  
  # profile is calculated from certain x location
  if isinstance(xLoc,int):
    if xLoc < 1024:
      maxHeight = np.max(np.mean(im[:,:950],axis=1))
    else:
      maxHeight = np.max(np.mean(im[:,1090:],axis=1))

    # averaging over several columns
    if averageFwhm>0:
      profile=np.median(im[:,xLoc-averageFwhm:xLoc+averageFwhm], axis=1)
    # only use a single column (could be affected by bad pixels)
    else:
      profile=im[:,xLoc]
    #TODO: lower limit for prominence of peaks is hard coded for now
    myProm = [50,maxHeight*3]
  else:
    if xLoc[0] < 1024:
      maxHeight = np.max(np.mean(im[:,:950],axis=1))
    else:
      maxHeight = np.max(np.mean(im[:,1090:],axis=1))
    profile=np.median(im[:,xLoc[0]:xLoc[1]], axis=1)
    #TODO: lower limit for prominence of peaks is hard coded for now
    myProm = [0,maxHeight*3]
#  print('max height for spatial fit is ', maxHeight)
  peaks = find_peaks(profile,width=peakWidth, prominence=myProm,height=[0,maxHeight*3])
  ind = peaks[0] 
  widths = peaks[1]['widths']
  prominences = peaks[1]['prominences']
  heights = peaks[1]['peak_heights']
  
  return ind, widths, prominences, heights, profile

def cnFindSpectralPeaks(im,peakWidth,spatialPeaks,spatialWidths,side):
  #TODO: implement subpixel width support
  assert(spatialWidths.dtype == np.int16()), "Spatial widths need to be integers"
  
  nrSpectra = len(spatialPeaks)
  profile = np.zeros((nrSpectra,1024))
  result = []
  if side == "left":
    for ii in range(nrSpectra):
      profile[ii,:] = np.median(im[spatialPeaks[ii]-spatialWidths[ii]:
        spatialPeaks[ii]+spatialWidths[ii],:1024],axis=0)
      peaks = find_peaks(profile[ii,:],width=peakWidth, prominence=[0.02,2000],height=[0,2000])
      ind = peaks[0] 
      widths = peaks[1]['widths']
      prominences = peaks[1]['prominences']
      heights = peaks[1]['peak_heights']
      result.append((ind,widths,prominences,heights))
  return result, profile

# total of 42 groups with an extra center pinhole
spacing = 14
center = np.arange(36,2048,47)
left = np.concatenate((center[0:21],center[22:]))+ spacing
right = np.concatenate((center[0:21],center[22:]))- spacing
y = np.arange(2048)

allY = np.concatenate((left,center,right))
profilePinhole = np.zeros(2048)
profileDiskSlit = np.zeros(2048)
profileDiskSlit[512:1536] = profileDiskSlit[512:1536]+1
kernel = gauss(np.arange(100)-50,1,0,20,0)
spatialDisk = np.convolve(profileDiskSlit,kernel,mode="same")/sum(kernel)
spatialFwhmPinhole = 3 # pixels

for ii in range(len(allY)):
  profilePinhole = profilePinhole + gauss(y,0.1,allY[ii],spatialFwhmPinhole/2.35,0)
  
siSpectrum = np.load("data/spectra/modulated-8-SiIX.npy")
tharSpectrum = np.load("data/spectra/SiIX-ThAr-spectrum.npy")
spectral = np.concatenate((tharSpectrum,tharSpectrum))
spectral = np.tile(spectral,(2048,1))+0.01

imPinhole = np.tile(profilePinhole,(2048,1)).T * spectral
imPinhole = rotate(imPinhole,0)
di = int((imPinhole.shape[0]-2048)/2)
imPinhole = imPinhole[di:di+2048,di:di+2048]
imPinhole = imPinhole/np.max(imPinhole)
#fig, ax=plt.subplots()
#ax.imshow(imPinhole,vmin=0,vmax=1)

imDisk = np.tile(spatialDisk, (2048,1)).T* spectral
imDisk = imDisk/np.max(imDisk)

im = imDisk
fig, ax=plt.subplots()
ax.imshow(im,vmin=0,vmax=1)

# For emmission spectra use average over almost all columns of beam
ind, widths, prominences, heights, profile= cnFindSpatialPeaks(im,[2,20],
                                                               xLoc=[10,930],
                                                               averageFwhm=None)
#fig, ax=plt.subplots()
#ax.plot(y,profile, y[ind], profile[ind], 'ro')

specResult, specProfile= cnFindSpectralPeaks(im,[2,20],ind,np.int16(widths),"left")

wo =0
specProfile = specProfile[wo,:]
centers = specResult[wo][0]
sigmas = specResult[wo][1]
proms = specResult[wo][2]
peakHeights = specResult[wo][3]
continuum = np.mean(peakHeights-proms) 

wl = np.arange(1024)
fit = np.zeros(wl.shape) 
for ii in range(len(centers)):
  fit = fit + gauss(wl,proms[ii],centers[ii],sigmas[ii]/(2.35),0)
fit = fit + continuum
#fig, ax=plt.subplots()
#ax.plot(wl,specProfile, wl, fit, 'r')

# some real data now
#real = fits.open("23marfocus_-40000.fits")[0].data.astype(np.float32)
#wind, widths, prominences, heights, profile= cnFindSpatialPeaks(-1*real,[1,10],
#                                                                xLoc=512,
#                                                                averageFwhm=50)
#print("# of peaks should be )",len(wind))
#fig, (ax1,ax2)=plt.subplots(2)
#ax1.plot(y,profile,y[wind],profile[wind],'ro')
#ax1.set_ylim([0,np.max(profile)])
#ax2.plot(wind,widths)
#ax2.set_ylim([1,4])
#
## do iteratively for every column
#my = [] 
#fig, ax=plt.subplots()
#ax.imshow(-1*real,vmin=0,vmax=3000,aspect='auto')
#for ii in np.arange(25,950):
#  ind, widths, prominences, heights, profile= cnFindSpatialPeaks(-1*real,[1,5],
#                                                                xLoc=int(ii),
#                                                                averageFwhm=25)
#  my.append((ind,widths))
##  print("# of peaks found(should be about 120)",len(ind))
#  ax.plot(np.ones(len(ind))*ii,ind,'r.',
#          np.ones(len(ind))*ii,ind-widths,'g.',
#          np.ones(len(ind))*ii,ind+widths,'m.',)
#   
#fig, ax=plt.subplots()                                                          
#for ii in np.arange(25,len(my),50):
#  ax.plot(my[ii][0],my[ii][1],label=str(ii))
#ax.legend()
plt.show()