#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
    
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *

def calWavecal(data,
               gratingPosition,
               dark,
               gain,
               badPixels,
               oldWavecal,
               simulateChange=False,
               threshold=None,
               mode=None,
               writeToFile=False,
               path=None,
               sequenceName=None,
               fileFormat='fits'):
  
  """
  Returns wavlength calibration table for the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    gratingPosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the grating stage for each ramp.
    dark : (#NDRs, 2048, 2048) ndarray, float32
        stored background dark master ramp from calibration store
    gain: (2048, 2048) ndarray, float32
        stored gain table from calibration store
    badPixels: (2048, 2048) ndarray, unit16
        stored bad pixel mask
    oldWavecal: (2048, 2048) ndarray, float32
        previous wavelength calibration from calibration store
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
    gainTable : (2048, 2048) ndarray, float32
        wavelength calibration table
    changeFlag: bolean
        indication whether wavecal has changed

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
  # TODO: helper function to find spatial bands for pinhole mask
  # TODO: helper function to find spectral lines in pinhole spectra
  # TODO: helper function for gaussian fit to spatial or spectral profile
  # TODO: remapping the beams would greatly help with detection of edges an fitting
  
  # Sixth determine if focus has changed
  # TODO: what is the threshold for a change
  
  if simulateChange:
    newWavecal = oldWavecal+10.0 # assuming user units simulate change by 10nm
    changeFlag = True
  else:
    newWavecal = oldWavecal
    changeFlag = False
  
  
  
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == 'fits':
      hdu = fits.PrimaryHDU(newWavecal)
      hdu.writeto(path+sequenceName+'.fits',overwrite=True)
    else:
      newWavecal.tofile(path+sequenceName+'.arr',sep="")
  
  return newWavecal, changeFlag


#a=cnH2rgRamps("data/gain/simGainVariation*",readMode="SLOW",subArray=None,verbose=True)
#data = np.squeeze(a.read(dtype=np.uint16))
#linData=np.zeros(data.shape)
#for i in range(3):
#  linData[i,:,:,:] = linearity(data[i,:,:,:], threshold=0, mode="SLOW")
#c= masterGain(data, threshold=0, mode="SLOW", writeToFile=True, path='data/gain/',
#                     sequenceName='gainTable',fileFormat='fits')

#fig, ax=plt.subplots()
#plt.imshow(data[0,5,:,:])
#
#fig, ax=plt.subplots()
#plt.imshow(data[0,5,:,:]*c)
#plt.show()