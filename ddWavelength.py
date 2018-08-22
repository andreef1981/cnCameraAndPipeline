#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *
#from helperFunctions import *
#from astropy.io import fits  # Reading/writing FITS data

def wavelength(data,
         waveCal,
         beamMapping):
  """
  Returns an array containing the wavelength for each pixel CryoNIRSP H2RG. 
  I keep the beamMapping information from ALIGN tasks here so we could
  use it to do the mapping of the two beams. This might either be done before
  or after demodulation.
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be mapped corrected.
    waveCal : (2048, 2048) ndarray, float32
        output from WAVECAL calibration task describing what the wavelength
        at each pixel is
    beamMapping : (2, 2048, 2048) ndarray, float32
        output from ALIGN3 calibration task
        information describing how left and right side of beam match. 

    Returns
    -------
    waveVector : (2048,) ndarray, float32
        vector containing the wavelength at each pixel

    Raises
    ------

    Notes
    -----
    Potentially we might want to remove curvature of lines in this and also do
    some beam mapping
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> 
   """
  
  # for now just use the center row
  waveVector = waveCal[1024,:]
  
  return waveVector

# sample code to run wavecal plugin
#a=cnH2rgRamps("data/gain/simGain*",readMode="SLOW",subArray=None,verbose=True)
#data = np.squeeze(a.read(dtype=np.uint16))
#
#waveCal = fits.open("data/wavecal/wavecal.fits")[0].data.astype(np.float32)
#b = wavelength(data,waveCal,1)
#
#plt.plot(b)