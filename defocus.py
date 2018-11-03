# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
#from datetime import datetime
from astropy.io import fits
from astropy.convolution import convolve
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
from numpy.linalg import linalg

siSpectrum = np.load("data/spectra/modulated-8-SiIX.npy")

test = np.squeeze(siSpectrum[0,:])
widths = np.array([20,15,10,5,1,6,11,16,21])
result = np.zeros((len(widths),len(test)))
for i in range(len(widths)):
  kernel = cnGauss(np.arange(-25,26),1,0,widths[i],0)
  result[i,:] = convolve(test,kernel,boundary='extend',normalize_kernel=True)
  
#%%
x = np.arange(2048)
fig, ax=plt.subplots(num=1)
ax.plot(x,test, x,result[4,:],'g',result[3,:],'r')
#%%
# find some peaks
ind = []
mw = []
prom = []
mi = np.zeros(len(widths),dtype=np.int16())
av = np.zeros(len(widths),dtype=np.float16())
for i in range(len(widths)):
  spectralPeaks = cnFindSpectralPeaks(result[i,:1024],
                                    [2,100],
                                    None,
                                    None,
                                    None,
                                    [0.05,1],
                                    [-1,0],
                                    invert=True)
  mi[i] = np.argmax(spectralPeaks[0][2])
  av[i] = np.mean(spectralPeaks[0][1])
  print(spectralPeaks[0][0][mi[i]])
  ind.append(spectralPeaks[0][0][mi[i]])
  mw.append(spectralPeaks[0][1][mi[i]])
  prom.append(spectralPeaks[0][2])
  
# for now use just the strongest peak
print(av)
#%%
# fig, ax=plt.subplots(num=1)
# ax.plot(x,-1*result[4,:],x[ind[4]],-1*result[4,ind[4]], 'ro')
# xx = np.array([-7.8,-9, -9.6, -9.9, -10.1, -11.3,-11.6,-12.2,-13.4])
xx = np.array([-9.4, -9.6, -9.8, -10., -10.2, -10.4, -10.6, -10.8, -11.])
par = np.polyfit(xx,np.float64(av),2)
fit = np.polyval(par,xx)
hr = np.arange(-11.4,-9,0.01)
fitter = np.polyval(par,hr)
wo = np.argmin(fitter)
print(hr[wo])
fig, ax=plt.subplots(num=1)
ax.plot(xx,av,xx,fit,'r', hr, fitter, 'g')



plt.show()