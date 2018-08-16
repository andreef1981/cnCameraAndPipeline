#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 12:07:17 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
15 August 2018:
    Removed instrument dark as it is contained in the background dark already.
    Andre
"""

import numpy as np
import matplotlib.pyplot as plt
from cnPipeline import *

def dark(data,,
         backgroundDark):
  
  #TODO: remove instrument dark subtraction and make sure background is linearity corrected
  """
  Returns the instrument and background dark subtracted ramp of the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#NDRs, 2048, 2048) ndarray, uint16
        3D data cube that needs to be dark corrected.
    backgroundDark : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that contains the background dark rampwith the exact same
        camera setup as the data cube. Has to be linearity corrected.

    Returns
    -------
    darkSubtracted : (#NDRs, 2048, 2048) ndarray, float32
        3D data cube that is dark corrected.

    Raises
    ------
    AssertationError
        If the shapes of the input arrays do not match.

    Notes
    -----
    
    CSS uses little endian (16 or 32 bits). BTD converts to big endian for 
    network transport. All simulated data should be little endian.
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> data = np.zeros((5,10,10),dtype='uint16')+6
    >>> instrumentDark = np.zeros((5,10,10),dtype='float32')+1
    >>> backgroundDark = np.zeros((5,10,10),dtype='float32')+2
    >>> darkSubtracted = dark(data,instrumentDark,backgroundDark)
   """
  
  #TODO: check input type, binary vs ndarray.

  
  #TODO: Should we implement a check for exposure time, camera mode or is the
  #      check for the right calibration data instrumentDark and backgroundDark
  #      done before calling plugin.
  
  
  # make sure the data cubes have the same dimensions
  assert(data.shape == backgroundDark.shape),\
  'input array dimensions do not match'
  
  #TODO: Do we need to check for data types and convert to float32?
#  if (instrumentDark.dtype is not 'float32' or backgroundDark.dtype is not 'float32'):
#    instrumentDark = np.float32(instrumentDark)
#    backgroundDark = np.float32(backgroundDark)
    
  
  # subtract the data
  darkSubtracted = data-backgroundDark
  
  return darkSubtracted

# sample code to run dark plugin
#a=cnH2rgRamps("data/instrumentDark/masterInstrumentDark*",readMode="SLOW",subArray=None,verbose=True)
#instrumentDark = np.squeeze(a.read(dtype=np.float32))
#
#b=cnH2rgRamps("data/backgroundDark/masterBackgroundDark*",readMode="SLOW",subArray=None,verbose=True)
#backgroundDark = np.squeeze(b.read(dtype=np.float32))
#c= cnH2rgRamps("data/backgroundDark/simBackgroundDark*",readMode="SLOW",subArray=None,verbose=True)
#data = c.read()[0,:,:,:]
#darkSubtracted = dark(data,instrumentDark,backgroundDark)
#
#fig, ax=plt.subplots()
#plt.imshow(data[9,:,:])
#fig, ax=plt.subplots()
#plt.imshow(darkSubtracted[9,:,:])
