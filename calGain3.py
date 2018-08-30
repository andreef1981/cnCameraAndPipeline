#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: Andre Fehlmann (afehlmann@nso.edu)

Revision history
----------------
28 August 2018:
  - rename the function to match others
  - add inputs: stage position, dark, badPixels, threshold, simulateChange
  - add ouputs: changeFlag
    
"""
import numpy as np
from astropy.io import fits
from cnPipeline import *
from ddLinearity import *

def calGain3(data,
             stagePosition,
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
  Returns the simple gain table for the spectrograph CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
    stagePosition : (#ramps,) ndarray, float32/uint64
        1D array that contains the positions of the lamp stage for each ramp.
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
  
  if len(data.shape) == 3:
    # only a single ramp is provided but matrix quadfit expects another axis
    data = np.expand_dims(data, axis=0)
    
  # matrixQuadfit needs float data to work properly
  data = np.float32(data)
    
  # in this case we use the NDR number as our 'time' axis
#  dataTime = np.arange(data.shape[1])
  
  # fit quadratic to data
  a,b,c = matrixQuadfit(data,threshold=threshold,mode=mode,ignoreRef=False)
  # Use numpy's multiply capability to to the multiplication along the right axis
#  nonLinear = np.multiply(np.expand_dims(linearityCoefficients[0,:,:],axis=0),
#                          dataTime[:,None,None]**2.)
  
  temp = np.mean(b,axis=0)
  gainTable = np.mean(temp)/temp
  
  # for test purposes lets write the files to fits
  if writeToFile:
    if fileFormat == 'fits':
      hdu = fits.PrimaryHDU(gainTable)
      hdu.writeto(path+sequenceName+'.fits',overwrite=True)
    else:
      gainTable.tofile(path+sequenceName+'.arr',sep="")
  
  return gainTable


#a=cnH2rgRamps("data/gain/simGainVariation*",readMode="SLOW",subArray=None,verbose=True)
#data = np.squeeze(a.read(dtype=np.uint16))
#c= masterGain3(data, threshold=0, mode="SLOW", writeToFile=True, path='data/gain/',
#                     sequenceName='gainTable',fileFormat='arr')
#
#fig, ax=plt.subplots()
#plt.imshow(data[0,5,:,:])
#
#fig, ax=plt.subplots()
#plt.imshow(data[0,5,:,:]*c)
#plt.show()