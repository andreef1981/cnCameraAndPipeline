#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 24 12:33:51 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt

beamMapping = np.ones((2048,2048),dtype="float32")
badPixels = np.ones((2048,2048),dtype="uint8")

data = beamMapping
file = "data/coronalObs-sensitivity/spBeamMapping"
#data = badPixels
#file = "data/coronalObs-sensitivity/ciBadPixels"
hdu = fits.PrimaryHDU(data)
hdu.writeto(file+'.{0:03d}'.format(1)+'.fits',overwrite=True)
data.tofile(file+'.{0:03d}'.format(1)+'.raw',sep="")