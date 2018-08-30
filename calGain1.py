#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:45:58 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

    
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *

def calGain1(data,
             stagePositionFM1A,
             stagePositionFM1B,
             dark,
             badPixels,
             oldGain,
             threshold=None,
             simulateChange=False):
  
  #TODO: are the stage positions in user or raw units
  #      (float32/uint64 for positions)
  #TODO: should we add a threshold in a property database to decide if gain
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the a simple gain table for the CI only.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best focus position.
    stagePositionFM1A : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the FM1A stage for each ramp.
    stagePositionFM1B : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the FM1B stage for each ramp.
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
    
    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> 
   """
  # First perform linearity correction
  
  # Second subtract linearized background
  
  # Third perform flat fielding
  
  # Fourth ignore/interpolate bad pixels
  
  # Fifth loop through sequences to find spatial and spectral focuses
  # TODO: helper function to find spatial bands for pinhole mask
  # TODO: helper function to find spectral lines in pinhole spectra
  # TODO: helper function for gaussian fit to spatial or spectral profile
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


