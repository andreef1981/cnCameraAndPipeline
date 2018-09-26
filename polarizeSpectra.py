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

def gauss(x,a,x0,sigma,c):
    return c+a*exp(-(x-x0)**2/(2*sigma**2))

#waveCal = fits.open("data/coronalObs-sensitivity/spWavecal.000.fits")[0].data.astype(np.float32)
#waveVector = waveCal[1024,:1024]

siVector = np.load("siWavelength.npy")
#heSpectrum = fits.open("data/spectra/he_spectrum_combined-32bitfloat.fits")[0].data.astype(np.float32)
#heVector = heSpectrum[1024,:1024]

siSpectrum = fits.open("data/spectra/siIX_spectrum_atmosphere.fits")[0].data.astype(np.float32)
siSpectrum = siSpectrum[1024,:1024]*1.06

print(siVector[512])
#wlaverage = np.median(subArr,axis=1)
#wlaverage = wlaverage+np.abs(np.max(wlaverage))
#  x = ar(range(np.size(wlaverage)))
#  y = ar(wlaverage)
wv0      = 3934.3
linewidth = 3e-4*wv0/(2*2.35) #FWHM to 2sigma to sigma
linestrength = 0.3
dop       = 0.10
azi       = np.deg2rad(315)
blos      = 10. ## Gauss
g=1.5
bfactor    = 4.67e-12*g*blos*wv0**2
print(bfactor)


#I = gauss(np.arange(len(siVector)), 1,512.,14.4,0)
#dI = np.multiply(I,((np.arange(len(siVector))-512.)*-1./(14.4**2)))

I = gauss(siVector, 1,wv0,linewidth,0)
dI = np.multiply(I,(siVector-wv0)*-1./(linewidth**2))
Q = dop * I * np.cos(azi)
U = dop * I * np.sin(azi)
V = bfactor*dI
print(np.max(dI),np.min(dI))
print(np.max(Q),np.min(Q))
print(np.max(V),np.min(V))
S = np.vstack((linestrength*I+siSpectrum,Q,U,V))

mod = np.load("8state_modmatrix.npy")
demod = np.load("8state_demodmatrix.npy")

obs = np.matmul(mod,S)
recov = np.matmul(demod,obs)
#comp = np.zeros((8,1024))
#for i in range(8):
#  for j in range(1024):
#    comp[i,j] = np.dot(mod[i,:],S[:,j])
#
#so = test-comp  
#print(np.min(so), np.max(so))
#popt,pcov = curve_fit(gauss,np.arange(len(heVector)),heVector,p0=[-0.052,519.,7.7,1.029])
#print('a = ',popt[0])
#print('center = ',popt[1])
#print('sigma = ',popt[2])
#print('FWHM = ',2.35*popt[2])

# figure 1
fig, (ax1,ax2)=plt.subplots(2,1)
fig.suptitle("emmission line and atmospheric spectrum")
ax1.plot(siVector, siSpectrum, siVector, I)
ax2.plot(siVector, siSpectrum+linestrength*I)

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)
fig.suptitle("pure emmission line Stokes vectors")
ax1.plot(siVector, I, siVector, dI)
ax2.plot(siVector, Q)
ax3.plot(siVector, U)
ax4.plot(siVector, V)


fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)
fig.suptitle("Stokes vectors with atmospheric contribution to I")
ax1.plot(siVector, siSpectrum+linestrength*I)
ax2.plot(siVector, Q)
ax3.plot(siVector, U)
ax4.plot(siVector, V)

fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(8,1)
fig.suptitle("observed (modulated) intensities")
ax1.plot(siVector, obs[0,:])
ax2.plot(siVector, obs[1,:])
ax3.plot(siVector, obs[2,:])
ax4.plot(siVector, obs[3,:])
ax5.plot(siVector, obs[4,:])
ax6.plot(siVector, obs[5,:])
ax7.plot(siVector, obs[6,:])
ax8.plot(siVector, obs[7,:])

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4,1)
fig.suptitle("observed (modulated) intensities")
ax1.plot(siVector, recov[0,:])
ax2.plot(siVector, recov[1,:])
ax3.plot(siVector, recov[2,:])
ax4.plot(siVector, recov[3,:])


im = np.concatenate((np.tile(siSpectrum+linestrength*I,(2048,1)),
                         np.tile(siSpectrum+linestrength*I,(2048,1))),
                        axis=1)
fig3, ax1 = plt.subplots()
plt.imshow(im)
fig3.tight_layout()
plt.show()