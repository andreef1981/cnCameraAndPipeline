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

def calFocus1(data,
              dark,
              gain,
              badPixels,
              oldFocus,
              simulateChange=False):
  
  #TODO: should we return the focus position in user or raw units
  #TODO: should we add a threshold in a property database to decide if focus
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the best focus position for the internal spectrograph focus SM5.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best focus position.
    dark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gain: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldFocus: float32/int64
        the prevously determined focus position
    simulateChange: bolean
        test flag

    Returns
    -------
    newFocus : float32/int64
        the new focus position for SM5
    changeFlag: bolean
        indication whether focus has changed

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
    newFocus = oldFocus+1.0 # assuming user units simulate change by 1mm
    changeFlag = True
  else:
    newFocus = oldFocus
    changeFlag = False
  
  return newFocus, changeFlag


