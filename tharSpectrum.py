#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 11:08:10 2018

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

def gauss(x,a,x0,sigma,c):
    return c+a*exp(-(x-x0)**2/(2*sigma**2))
  
def wnum2ang(wavenum, help=False):
 
  if help is True:
    print('Convert input wavenumber in cm^-1 (vacuum) to air Angstroms')
    print('angstrom=wnum2ang(wavenum)')
    print('wavenum - floating point wavenumber in Kaysers')
    print('(cm^-1 in vacuum)')
    print('angstrom - wavelength in Angstroms computed for air')
    angstrom=0.
    return angstrom
	

  #Compute index of refraction from NBS formulation
  a=1.+6432.8e-8
  b=2949810.0
  c=146.0e8
  d=25540.0
  e=41.0e8
  n=a+(b/(c-wavenum**2.))+(d/(e-wavenum**2.))
  
  angstrom=1./(n*wavenum)*1.e8
  return angstrom

def tharSpectrum(cw,
                 wavelengthOnArray,
                 bandpass):
#  3934 nm
  wavenumber = np.array([2530.4232, 2531.7493, 2542.4886, 2542.8614, 2549.9993, 2552.3161])
  waveWidth = np.array([0.014, 0.021, 0.014, 0.025, 0.016, 0.025])#fwhm = 2.355*sigma
  intensity = np.array([2,20,9,53,3,32])
  names = ["Th I", "Ar I", "Th II", "Ar I", "Th I", "Ar I"]
  
  # 1083.0 nm
  #Engleman 2003
#  wavenumber = np.array([9136.3890, 9136,9082, 9151.2409, 9160.3878, 9162.2777,
#                         9170.7994, 9179.9819])
#  waveWidth = np.array([0.020,0.019,0.018,0.021,0.022,0.026,0.028])
#  intensity = np.array([4215,296,111,98,80,410,464])
#  names = ["Th II", "Th II", "Th II", "Th I", "Th I", "Th I", "Th I"]
  # Hinkle 2001
#  scale = 1.#np.mean(np.array([80/0.0066,98/0.0191,410/0.0326,464/0.0243]))
#  wavenumber = np.array([9152.115,9160.388,9162.278,
#                         9170.799,9173.53,9175.288,9176.116,9178.231,9179.982,
#                         9183.69,9185.371,9186.487,9187.864,9189.64,9197.142,
#                         9198.903,9199.360,9199.695,9201.688,9203.462,9204.055,
#                         9204.668,9207.824,9209.04,9211.130,9217.941,9218.820,
#                         9221.51,9222.53,9224.63,9226.583,9227.152,9227.527,
#                         9229.46,9231.543,9235.95,9236.23,9237.31,9239.545,
#                         9241.423,9242.606,9245.257,9245.683,9246.209,9250.433,
#                         9250.707,9252.848,9253.365,9254.553,9256.682,9258.28,
#                         9258.75,9260.11,9261.815,9268.425,9268.815,9269.478,
#                         9274.675,9276.494,9279.606,9282.150])
#  waveWidth = np.ones(wavenumber.shape)*0.6#/(2.35)
#  intensity = scale*np.array([0.0610,0.0191,0.0066,
#                        0.0326,0.0015,0.0013,0.0023,0.2303,0.0243,0.0101,0.0024,
#                        0.0077,1.5511,0.0010,0.0415,0.0063,0.0156,0.0053,0.0064,
#                        0.0230,0.0017,0.1426,0.0065,0.0017,0.0196,0.0072,0.0250,
#                        0.0016,0.0026,0.0013,0.0039,0.0039,0.0111,0.0017,0.0122,
#                        0.0024,0.0026,0.0097,0.0786,0.0078,0.0256,0.4545,0.1455,
#                        0.0626,0.0192,0.0472,0.0034,0.1019,0.0038,0.1442,0.0013,
#                        0.0023,0.0055,0.0111,0.0121,0.0557,0.0079,0.0080,0.0526,
#                        0.4566,0.0736])
#  names = ["Ar II","Th I","Th I",
#           "Th I","Ar II","Ar I","Th I","Ar I","Th I","Ar I","Th I","Th I",
#           "Ar I","Undef","Th I","Ar II","Ar II","Th I","Th I","Th I","Th II",
#           "Ar I","Th I","Undef","Th I","Ar I","Th I","Undef","Ar II","Undef",
#           "Th blend","Th II","Th I","Ar I","Ar II","Ar I","Undef","Ar I",
#           "Th II","Ar II","Ar II","Th I","Ar II","Th I","Th I","Ar I","Th I",
#           "Th I","Th I","Th I","Undef","Undef","Ar I","Th II","Th I","Th I",
#           "Ar II","Th I","Th I","Ar I","Ar I"]

  wavelength = wnum2ang(wavenumber)/10.
  
  fwhmWavelength = waveWidth/wavenumber*wavelength
  
  artRes = 0.001
  # dispersion in nm / pixel
  dispersion = (np.max(wavelengthOnArray)-np.min(wavelengthOnArray))/1024.
  print(dispersion)
  slitProfile = np.int(np.ceil(6*dispersion/artRes))
  print(slitProfile)
  x = np.arange(-200,200,1)
  kernel = gauss(x,1,0,slitProfile/(2.35),0)
#  kernel = np.float32(np.zeros(slitProfile+2))
#  kernel[1:-1]=1
  wl = np.arange(cw-bandpass,cw+bandpass,artRes)
  spec = np.zeros(wl.shape)
  
  for ii in range(len(wavenumber)):
    spec = spec + gauss(wl,intensity[ii],wavelength[ii],fwhmWavelength[ii]/(2.35),0)
  spec = spec/np.max(spec)
  conv = np.convolve(spec,kernel,mode="same")/sum(kernel)
  # interpolate to 
  newx = np.arange(np.min(wavelengthOnArray), np.max(wavelengthOnArray), dispersion)
  allx = np.arange(np.min(wl), np.max(wl), dispersion)
  final = np.interp(newx,wl,conv)
  allfinal = np.interp(allx,wl,conv)
  return wl, spec, conv, newx, final,allx,allfinal



cw = 3934.3
wavelengthOnArray= np.array([3925.1, 3942.9])
#cw = 1083.0
#wavelengthOnArray= np.array([1080.8, 1085.2])
bandpass = 20

wavelength, spec,conv, newx,final,allx,allfinal = tharSpectrum(cw,wavelengthOnArray,bandpass)
#np.save("data/spectra/SiIX-ThAr-wavelength.npy", newx)
#np.save("data/spectra/SiIX-ThAr-spectrum.npy", final)

refSpectrum = np.concatenate((allx[None,],allfinal[None,]),axis=0)
# np.save("data/spectra/SiIX-ThAr-reference-spectrum-175um.npy", refSpectrum)
fig, ax=plt.subplots()
ax.plot(wavelength,spec,wavelength,conv,'r', newx,final,'g')

allfinal = 30*allfinal

test=np.correlate(allfinal,final,mode="same")
print(np.max(test))
fig, ax=plt.subplots()
ax.plot(test)
ind = np.where(test==np.max(test))
fig, ax=plt.subplots()
ax.plot(allx[ind[0][0]-512:ind[0][0]+512],allfinal[ind[0][0]-512:ind[0][0]+512],'o',
        newx,final)

