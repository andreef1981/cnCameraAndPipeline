#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

"""
import numpy as np

def calAlign1(data,
              stagePositionFM1A,
              stagePositionFM1B,
              camera,
              dark,
              gain,
              badPixels,
              oldAlignFM1A,
              oldAlignFM1B,
              oldStepSizeFM1A,
              oldStepSizeFM1B,
              oldRotationFM1A,
              oldRotationFM1B,
              threshold=None,
              simulateChange=np.zeros((2,3),dtype=bool)):
  
  #???: should we return the align position in user or raw units
  #      (float32/uint64 for positions)
  #???: should we add a threshold in a property database to decide if align
  #      position has changed?
  #???: Where do we keep the bad pixel mask stored?
  #???: where is delta calculation to base done?
  #???: does the dark need to be a "GOS thermal dark"?
  
  """
  Returns the best center position of the beam steering mirror FM1A and FM1B 
  for a given slit mask, pickoff state. 
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best alignment.
        (raster scan)
    stagePositionFM1A : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the FM1A stage for each ramp.
    stagePositionFM1B : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the FM1B stage for each ramp.
    camera : string
        string describing which instrument is used "SP" or "CI"
    dark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gain: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldAlignFM1A: float32/int64
        the previously determined best alignment for FM1A
    oldAlignFM1B: float32/int64
        the previously determined best alignment for FM1B
    oldStepSizeFM1A: float32
        the previously determined step size for FM1A
    oldStepSizeFM1B: float32
        the previously determined step size for FM1B
    oldRotationsFM1A: float32
        the previously determined rotation for FM1A
    oldRotationFM1B: float32
        the previously determined rotation for FM1B
    threshold : (2,3) ndarray, float32/uint64
        change beyond which change flag is set
    simulateChange: (2,3) ndarray, boolean
        boolean array to simulate change in alginments

    Returns
    -------
    newAlignFM1A : float32/int64
        the new best alignment for FM1A
    newAlignFM1B : float32/int64
        the new best alignment for FM1B
    newStepSizeFM1A: float32
        the previously determined step size for FM1A
    newStepSizeFM1B: float32
        the previously determined step size for FM1B
    newRotationsFM1A: float32
        the previously determined rotation for FM1A
    newRotationFM1B: float32
        the previously determined rotation for FM1B
    changeFlag: (2,3) ndarray, boolean
        boolean array to indicat changes in FM1A or FM1B

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
   # initialize the change flag to all False
  changeFlag=np.zeros((2,3),dtype=bool)
   
  # First perform linearity correction
  
  # Second subtract linearized background
  
  # Third perform flat fielding
  
  # Fourth ignore/interpolate bad pixels
  
  # Fifth decide which instrument is used and loop through sequences to find
  # best alignment
  # TODO: helper function to find brightest spectra for SP (what is the GOS
  #       pinhole image size on the SP slit plane?) that is centered on slit
  #       image on array (data is raster scan)
  # TODO: for SP determine rotation of slit relative to steering mirror axes.
  # TODO: for CI; helper function to find FM1 positions that center
  #       GOS pinhole on array 
  # TODO: for CI determine rotation of array relative to steering mirror axes.
  # TODO: for both calibrate step size to pixels
  
  # Sixth determine if alignment has changed
  # TODO: what is the threshold for a change
  
  if simulateChange[0,0]:
    newAlignFM1A = oldAlignFM1A+0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[0,0] = True
  if simulateChange[1,0]:
    newAlignFM1B = oldAlignFM1B-0.05 # assuming user units simulate change by 0.05 degrees
    changeFlag[1,0] = True
  if simulateChange[0,1]:
    newStepSizeFM1A = oldStepSizeFM1A+3.0 # simulate change by 3 mu per pixel
    changeFlag[0,1] = True
  if simulateChange[1,1]:
    newStepSizeFM1B = oldStepSizeFM1B-3.0 # simulate change by 3 mu per pixel
    changeFlag[1,1] = True
  if simulateChange[0,2]:
    newRotationFM1A = oldRotationFM1A+0.8 # simulate change by 0.8 degrees
    changeFlag[0,2] = True
  if simulateChange[1,2]:
    newRotationFM1B = oldRotationFM1B-0.8 # simulate change by 0.8 degrees
    changeFlag[1,2] = True
  
  return newAlignFM1A, newAlignFM1B, newStepSizeFM1A, newStepSizeFM1B,\
         newRotationFM1A, newRotationFM1B, changeFlag

