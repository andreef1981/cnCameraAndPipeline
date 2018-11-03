#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 08:52:45 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
from helperFunctions import *
from datetime import datetime
from astropy.io import fits
from matplotlib.colors import LogNorm, PowerNorm
import sunpy.cm
import matplotlib.image as mpimg
from scipy import interpolate

# img = mpimg.imread('latest_1024_0193.jpg')

# xdim = img.shape[0]
# ydim = img.shape[1]
# x = np.arange(xdim)
# y = np.arange(ydim)
# newx = np.arange(0,xdim,2)
# g0 = interpolate.interp2d(x, y, img[:,:,0], kind='cubic')
# g1 = interpolate.interp2d(x, y, img[:,:,1], kind='cubic')
# g2 = interpolate.interp2d(x, y, img[:,:,2], kind='cubic')
# newIm = np.zeros((int(xdim/2),int(ydim/2),3))
# newIm[:,:,0] =g0(newx,newx)
# newIm[:,:,1] =g1(newx,newx)
# newIm[:,:,2] =g2(newx,newx)
# newIm = np.int16(newIm)
#mpimg.imsave('test.png', newIm, vmin=None,
#                        vmax=None, cmap=None,
#                        format=None, origin=None, dpi=300)
mpimg.thumbnail('latest_1024_0193.jpg', 'thumb.png',
                scale=0.5, interpolation='bilinear', preview=False)

# fig, ax=plt.subplots(num=1)
# cax=ax.imshow(newIm[:,:,0])
# fig.colorbar(cax)

#%%
hdul = fits.open("f0171.fits")#[0].data.astype(np.float32)

hdul.verify('silentfix')
# im, hdr = fits.getdata("f0171.fits", 1, header=True)
# hdul.info()
im = hdul[1].data
# print(np.min(im), np.max(im))
# hdr = fits.getheader("f0171.fits",1)
# im = fits.getdata('f0171.fits',1)
h = hdul[1].header
# print(repr(h))
# print(h['BITPIX'])
# hdr = h
im = np.array(im)
x=np.arange(4096)
scale = 8
newx=np.arange(0,4096,scale)
ycen = h['CRPIX1']
ysca = h['CDELT1'] #arcsec/pixel
xcen = h['CRPIX2']
xsca = h['CDELT2'] #arcsec/pixel
dsunRef = h['DSUN_REF'] # in m
dsunObs = h['DSUN_OBS'] # in m
rsun = h['R_SUN'] # in pixels (observed?)
rsunRef = h['RSUN_REF'] # in m
rsunObs = h['RSUN_OBS'] # arcsec
newxcen = np.int16(np.round(xcen/scale))
newycen = np.int16(np.round(ycen/scale))
rsunScale = rsunObs/959.5 # arcsec
rsunPixel = 959.5/(xsca*scale)
print(rsunPixel)


f = interpolate.interp2d(x, x, im, kind='cubic')
newIm = f(newx,newx)

finalIm = np.zeros((697,697))
finalIm[348-256:348+256,348-256:348+256]=newIm
finalIm = finalIm*10
finalIm = np.where(finalIm<1.,1,finalIm)
finalIm = np.int16(finalIm)
# finalIm = np.where(finalIm>400,4*finalIm,finalIm)
circle1 = plt.Circle((348, 348), rsunPixel, color='r',fill=False)
circle3 = plt.Circle((348, 348), rsunScale*rsunPixel, color='b',fill=False)
circle2 = plt.Circle((348, 348), 1.5*rsunPixel, color='g',fill=False)

hdul.close()
fig, ax=plt.subplots(num=1)
cax=ax.imshow(im)
ax.plot(xcen,ycen,'ro')
fig.colorbar(cax)

fig, ax=plt.subplots(num=2)
cax=ax.imshow(newIm)
ax.plot(newxcen,newycen,'ro')
fig.colorbar(cax)
#%%
fig, ax=plt.subplots(num=3)
cax=ax.imshow(finalIm, cmap=plt.get_cmap('afmhot'),norm=LogNorm(1000,
              vmax=np.max(finalIm)),interpolation='None')
# cax=ax.imshow(finalIm, cmap=plt.get_cmap('sdoaia171'),norm=LogNorm(1000,
#               vmax=np.max(finalIm)))
ax.plot(348,348,'ro')
ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
fig.colorbar(cax)

mpimg.imsave('test.png', finalIm, vmin=100,
              vmax=np.max(finalIm), cmap=plt.get_cmap('sdoaia171'),
                        format=None, origin=None, dpi=300)
plt.show()

