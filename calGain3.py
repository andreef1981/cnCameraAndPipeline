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

def masterGain3(data,
               threshold=None,
               mode=None,
               writeToFile=False,
               path=None,
               sequenceName=None,
               fileFormat='fits'):
  
  """
  Returns the gain table for the CryoNIRSP H2RG.
  
   Parameters
    ----------
    data : (#ramps, #NDRs, 2048, 2048) ndarray, float32
        4D data cube that will be averaged.
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
        gain table for multiplicative gain correction

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
  dataTime = np.arange(data.shape[1])
  
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