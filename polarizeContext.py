#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 11:37:03 2018

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
from matplotlib.patches import Circle


def gauss(x,a,x0,sigma,c):
    return c+a*exp(-(x-x0)**2/(2*sigma**2))

#waveCal = fits.open("data/coronalObs-sensitivity/spWavecal.000.fits")[0].data.astype(np.float32)
#waveVector = waveCal[1024,:1024]

def circleMask(index,radius,array):
  a,b = index
  ri,ro = radius
  nx,ny = array.shape
  y,x = np.ogrid[-a:nx-a,-b:ny-b]
  mask = (x*x + y*y <= ro*ro) & (x*x + y*y >= ri*ri) 
  array[mask]=1
  return array

def cross45Mask(distance,array):
  nx,ny = array.shape
  result = np.zeros([nx,ny])
  for i in range(distance):
    vec1 = np.arange(nx-i)
    temp1 = np.diagflat(vec1,k=i)
    temp2 = np.diagflat(vec1,k=-i)
    result = result + temp1 + temp2
  result = result + np.rot90(result)
  result[result>0] = 1
  return result

def crossMask(distance,array):
  
  
  nx,ny = array.shape
  result = np.zeros([nx,ny])
  result[:,int(nx/2)-1-distance-1:int(nx/2)+distance-1]=1
  result[int(ny/2)-1-distance-1:int(ny/2)+distance-1,:]=1
  result[result>0] = 1
  return result

circle = circleMask([1024,1024],[900,1000], np.zeros([2048,2048]))
cross = crossMask(50,np.zeros([2048,2048]))
cross45 = cross45Mask(100,np.zeros([2048,2048]))

names = ["latest_4096_0193-upper-left.png",
         "latest_4096_0193-upper-center.png",
         "latest_4096_0193-upper-right.png",
         "latest_4096_0193-center-left.png",
         "latest_4096_0193-center-center.png",
         "latest_4096_0193-center-right.png",
         "latest_4096_0193-lower-left.png",
         "latest_4096_0193-lower-center.png",
         "latest_4096_0193-lower-right.png"]
allObs = np.zeros([72,2048,2048])

for i in range(len(names)):

  im = plt.imread(names[i])
  im = im/np.max(im)
  print(im.shape)

# circle
# fig, ax = plt.subplots(num=1)
# ax.imshow(im/np.max(im))

# fig, ax = plt.subplots(num=2)
# ax.imshow(circle)

# fig, ax = plt.subplots(num=3)
# ax.imshow(cross)
#%%
  wv0      = 3934.3
  azi       = np.deg2rad(30)
  dop=0.1
  I = im
  Q = dop * I * np.cos(azi)*cross
  U = dop * I * np.sin(azi)*cross45
  V = 0.1*dop*I*circle

# fig, ax = plt.subplots(num=1)
# ax.imshow(I)

# fig, ax = plt.subplots(num=2)
# ax.imshow(Q)

# fig, ax = plt.subplots(num=3)
# ax.imshow(U)

# fig, ax = plt.subplots(num=4)
# ax.imshow(V)

  # now each pixel is going to be a linear combination of the 4 stokes signals
  fI = I.flatten()
  fQ = Q.flatten()
  fU = U.flatten()
  fV = V.flatten()
  S = np.vstack((fI,fQ,fU,fV))
  
  mod = np.load("8state_modmatrix.npy")
  demod = np.load("8state_demodmatrix.npy")
  
  obs = np.matmul(mod,S)
  recov = np.matmul(demod,obs)
  
  obsIm = np.reshape(obs,[-1,2048,2048])
  recovIm = np.reshape(recov,[-1,2048,2048])
  allObs[i*8:i*8+8,:,:] = obsIm
#%%
# fig, ax = plt.subplots(num=1)
# ax.imshow(obsIm[0,:,:])

# fig, ax = plt.subplots(num=2)
# ax.imshow(obsIm[1,:,:])

# fig, ax = plt.subplots(num=3)
# ax.imshow(obsIm[2,:,:])

# fig, ax = plt.subplots(num=4)
# ax.imshow(obsIm[3,:,:])

# fig, ax = plt.subplots(num=5)
# ax.imshow(obsIm[4,:,:])

# fig, ax = plt.subplots(num=6)
# ax.imshow(obsIm[5,:,:])

# fig, ax = plt.subplots(num=7)
# ax.imshow(obsIm[6,:,:])

# fig, ax = plt.subplots(num=8)
# ax.imshow(obsIm[7,:,:])

# fig, ax = plt.subplots(num=1)
# ax.imshow(recovIm[0,:,:])

# fig, ax = plt.subplots(num=2)
# ax.imshow(recovIm[1,:,:])

# fig, ax = plt.subplots(num=3)
# ax.imshow(recovIm[2,:,:])

# fig, ax = plt.subplots(num=4)
# ax.imshow(recovIm[3,:,:])



# plt.show()

#%%

np.save("data/spectra/modulated-8-9pos-contextImager.npy", allObs)
np.save("data/spectra/demod-8-SiX-CI.npy",demod)

