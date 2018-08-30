#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)
"""

import numpy as np
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *

def calGain4(data,
             gratingPosition,
             dark,
             badPixels,
             oldGain,
             threshold=None,
             simulateChange=False,
             mode=None,
             writeToFile=False,
             path=None,
             sequenceName=None,
             fileFormat='fits'):
  
  #TODO: are the stage positions in user or raw units
  #      (float32/uint64 for positions)
  #TODO: should we add a threshold in a property database to decide if gain
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the dither gain table for the spectrograph CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    gratingPosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the grating stage for each ramp.
    dark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldGain: (2048, 2048) ndarray, float32
        the prevously determined gain table
    threshold : float32/uint64
        change beyond which change flag is set
    simulateChange: bolean
        test flag
    writeToFile : bolean, optional, default=False
        writing to fits files for testing
    path: string
        directory to which optinal fits files will be written to
    sequenceName: string
        name of sequence to which optinal fits files will be written to
    fileFormat: string ('fits' | 'arr')
        writing fits files or little endian binary arrays

    Returns
    -------
    newGain: (2048, 2048) ndarray, float32
        new gain table for multiplicative gain correction
    changeFlag: bolean
        indication whether gain has changed

    Raises
    ------
    
    Notes
    -----
    
    CSS uses little endian (16 or 32 bits). BTD converts to big endian for 
    network transport. All simulated data should be little endian.
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> 
   """
  
  
  #TODO: reference pixel correction should be done prior to averaging. Need to check when and if to do it.
  
  # First perform linearity correction
  
  # Second subtract linearized background
  
  # Third perform flat fielding
  
  # Fourth ignore/interpolate bad pixels
  
  # Fifth loop through sequences to find spatial and spectral focuses
  # TODO: helper function to perform dither gain calibration
  # TODO: remapping the beams would greatly help with detection of edges an fitting
  
  # Sixth determine if focus has changed
  # TODO: what is the threshold for a change
  
  if simulateChange:
    newGain = oldGain*1.1 # simulate change by 10%
    changeFlag = True
  else:
    newGain = oldGain
    changeFlag = False
  
  return newGain, changeFlag
