#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 15:36:17 2018

@author: andreef
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *

import matplotlib.pyplot as plt

#a=np.arange(10)
#b = np.expand_dims(a - 3.2*a +1*a**2.,axis=1)
#c = np.repeat(b,1,axis=1)
#print(c[:,0])
#p= np.polyfit(a, c, 2)
#print(p)
#fit = np.polyval(p[:,0],a)
#fig, ax=plt.subplots()
#ax.plot(a,c[:,0],'o')
#ax.plot(a,fit,'r')
#plt.show()

#time = np.float16(np.arange(4))
#test = np.arange(100).reshape((5,5,4))
#test[0,0,:]=[1,2,3,4]
#new = test.reshape(-1, test.shape[-1])
#p = np.polyfit(time, new.T, 2)
#
#fit = np.polynomial.polynomial.polyval(time,np.flipud(p), tensor=True)

# reading byte arrays
#da = np.fromfile('data/thermalDark/masterThermalDark-00000-00009.arr',dtype=np.float32,count=-1,sep="")
#s = da.astype(bytes)
#st = s.astype(np.float32)
#da = np.memmap('data/thermalDark/masterThermalDark-00000-00009.arr',
#               dtype='float32', mode='r', shape=(2048,2048))
#fig, ax=plt.subplots()
#
#np.bytes_()
#plt.imshow(np.reshape(da,(2048,2048)))
##plt.imshow(da)
#plt.show()

#varVector = np.random.standard_normal(size=(10))
#varVector = np.random.normal(loc=1.0, scale=0.01, size=10)
#data = np.ones((10,10))
#res = np.multiply(data,varVector[:,None])

# He dispersion is 4.23e-3 nm per pixel
a = (np.arange(1024)-512)*4.23e-3 + 1083.
waveVector = np.concatenate((a,a))
waveMatrix = np.float32(np.repeat(waveVector[None,:],2048, axis=0))
hdu = fits.PrimaryHDU(waveMatrix)
hdu.writeto('data/wavecal/wavecal.fits',overwrite=True)

fig, ax=plt.subplots()
ax.plot(np.arange(2048), waveVector)