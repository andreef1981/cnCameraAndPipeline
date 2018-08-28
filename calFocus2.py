#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 09:46:37 2018

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------

"""
def calFocus2(data,
              stagePosition,
              instrument,
              dark,
              gain,
              badPixels,
              oldFocus,
              simulateChange=False):
  
  #TODO: should we return the focus position in user or raw units
  #      (float32/uint64 for positions)
  #TODO: should we add a threshold in a property database to decide if focus
  #      has changed?
  #TODO: Where do we keep the bad pixel mask stored?
  
  
  #TODO: where is delta calculation to base done?
  
  """
  Returns the best focus position for the external warm focus FM2.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, uint16
        4D data cube that will be used to determine the best focus position.
    stagePosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the focus stage for each ramp.
    instrument : string
        name string of instrument that will be focused
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
        the new focus position for FM2
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
  # TODO: for SP helper function to find spatial bands for GOS pinhole mask
  # TODO: for SP helper function to find spectral lines in GOS pinhole spectra
  # TODO: remapping the beams would greatly help with detection of edges an fitting
  # TODO: for CI helper function to find minimal spot size
  
  # Sixth determine if focus has changed
  # TODO: what is the threshold for a change
  
  if simulateChange:
    newFocus = oldFocus+1.0 # assuming user units simulate change by 1mm
    changeFlag = True
  else:
    newFocus = oldFocus
    changeFlag = False
  
  return newFocus, changeFlag

