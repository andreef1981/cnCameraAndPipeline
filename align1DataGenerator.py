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

myx = np.arange(-10, 10, 0.0001)

# spectrograph 0 order
xscale = 8.3333 #pixel/arcsec
yscale = 8.3333 #pixel/arcsec
xsize = 1024
ysize = 2048
width = 5.
givenX = np.arange(-19, 20, 2)
givenY = np.arange(-110, 120, 20)
x = np.int16(givenX*xscale + xsize/2 - 5)
y = np.int16(givenY*yscale + ysize/2)
# create all images
spZero = np.zeros([ysize,xsize,x.size,y.size],dtype=np.float16)
for j in range(y.size):
  for i in range(np.size(x)):
    print(i,j)
    # adding j tilts one axis
    spZero[:,:,i,j] = cnGauss2d(np.arange(xsize), np.arange(ysize), x[i]+j,
      y[j], width,width)
    
# context imager
xscale = 20 #pixel/arcsec
yscale = 20 #pixel/arcsec
xsize = 2048
ysize = 2048
width = 10.
givenX = np.arange(-45, 45, 5)
givenY = np.arange(-45, 45, 10)
x = np.int16(givenX*xscale + xsize/2 - 5)
y = np.int16(givenY*yscale + ysize/2)
# create all images
ciIm = np.zeros([ysize,xsize,x.size,y.size],dtype=np.float16)
for j in range(y.size):
  for i in range(np.size(x)):
    print(i,j)
    # adding j tilts one axis
    ciIm[:,:,i,j] = cnGauss2d(np.arange(xsize), np.arange(ysize), x[i]+j,
      y[j], width,width)
    
#%%    
fig, ax=plt.subplots(num=2)
ax.imshow(np.float32(ciIm[:,:,8,4]))

#%%
# base configuration for spectrograph in zero order similar to ci by finding center of spot
# first average in x and y direction to find center 
# givenX = np.arange(-10, 10, 1)+1/3.
# givenY = np.arange(-80, 80, 10)
# go through all images and find peak locations
xloc = np.zeros([x.size,y.size])
yloc = np.zeros([x.size,y.size])
xCenterPixel = np.zeros(y.size)
xStep = np.zeros(y.size)
for j in range(y.size):
  for i in range(np.size(x)):
    # print(i,j)
    xInd, xWidths, xProminences, xHeights, xProfile = cnFind2DSpatialPeak(im[:,:,i,j],
                                                                          0,
                                                                          [6,30],
                                                                          [0,ysize-1],
                                                                          None,
                                                                          [0.00001,10.],
                                                                          [0.0001,10.],
                                                                          invert=False)

    yInd, yWidths, yProminences, yHeights, yProfile = cnFind2DSpatialPeak(im[:,:,i,j],
                                                                          1,
                                                                          [6,30],
                                                                          [0,xsize-1],
                                                                          None,
                                                                          [0.00001,10.],
                                                                          [0.0001,10.],
                                                                          invert=False)
    
    # raise warning if one is bigger than size==1
    assert(xInd.size==1),"More than one peak in x direction found"
    assert(yInd.size==1),"More than one peak in y direction found"
    xloc[i,j] = xInd
    yloc[i,j] = yInd
  # do linear fit to find x scale and center at each y
  # p[0] is in arcsec/pixel
  p = np.polyfit(xloc[:,j],givenX,1)
  xCenterPixel[j] = -1*p[1]/p[0]
  xStep[j] = 1/p[0]

yCenterPixel = np.zeros(x.size)
yStep = np.zeros(x.size)
for i in range(np.size(x)):
  # do linear fit to find y scale and center at each x
  # p[0] is in arcsec/pixel
  p = np.polyfit(yloc[i,:],givenY,1)
  yCenterPixel[i] = -1*p[1]/p[0]
  yStep[i] = 1/p[0]
#%%
# need to find the arcsec value that puts spot on x=100 for y=100
xCenterOfY = np.zeros(y.size)
xSlopes = np.zeros(y.size)
for i in range(np.size(y)):
  p = np.polyfit(xCenterPixel,yloc[i,:],1)
  xCenterOfY[i] = (ysize/2-p[1])/p[0]
  xSlopes[i] = p[0]
xShiftArcsec = (np.mean(xCenterOfY)-xsize/2)/np.mean(xStep)
xAngle = np.rad2deg(np.arctan(xSlopes))

yCenterOfX = np.zeros(y.size)
ySlopes = np.zeros(x.size)
for i in range(np.size(y)):
  
  p = np.polyfit(xloc[:,i],yCenterPixel,1)
  yCenterOfX[i] = (xsize/2*p[0]) + p[1]
  ySlopes[i] = p[0]
yShiftArcsec = (np.mean(yCenterOfX)-ysize/2)/np.mean(yStep)
yAngle = np.rad2deg(np.arctan(ySlopes))
# print(yShiftArcsec)  
    

plt.show()

