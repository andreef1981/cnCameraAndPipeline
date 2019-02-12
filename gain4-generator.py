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
from logStepAlign3 import logStepAlign3

siVector = np.load("siWavelength.npy")
siSpectrum = fits.open("data/spectra/siIX_spectrum_atmosphere.fits")[0].data.astype(np.float32)
# siSpectrum = siSpectrum[1024,:1024]*1.06

#%%
test = np.squeeze(siSpectrum[0,:1024])
moves = np.arange(10)*10 - 40
result = np.zeros((len(moves),len(test)))
for i in range(moves.size):
  result[i,:] = np.roll(test,moves[i])

result = np.concatenate((result,1.06*result),axis=1)
np.save("data/spectra/rolled-10-SiIX.npy", result) 

#%%
fig, ax=plt.subplots(num=1)
ax.plot(result[4,:])



plt.show()